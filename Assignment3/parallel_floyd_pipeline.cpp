#include <iostream>
#include <iomanip>
#include <limits>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <algorithm>
#include "mpi.h"



// number of vertices
#define N 6

using namespace std;

double** alloc_2d_double(int rows, int cols);

void free_2d_double(double ** G, int rows, int cols);

void generate_graph(double ** G, int num_of_vertices);

void copy_2d_array(double ** destination, double ** source, int num);

void print_2d_array(double ** D, int rows, int cols);

void print_1d_array(double * a, int len);

int main(int argc, char* argv[]) {

    int id;   // process Id

    int p;    // number of processors

    int master_process = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    MPI_Status status;

    int d = sqrt(p);            // dimenstion, how many blocks in each dimension
    int bsz = N / sqrt(p);      // N/sqrt(p), block size in each process

    //****************** Distribute data into processes *****************************
    double ** block = alloc_2d_double(bsz, bsz);

    if (id == 0) {
        // root process generate Graph represented by 2d array
        double ** G = alloc_2d_double(N, N);
        generate_graph(G, N);
        print_2d_array(G, N, N); 

        // root process distribute data into other processes
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                double ** tmp_block = alloc_2d_double(bsz, bsz);

                for (int k = 0; k < bsz; k++) {
                    for (int l = 0; l < bsz; l++) {
                        tmp_block[k][l] = G[k + i * bsz][l + j * bsz];
                    }
                }

                int des_pid = i * d + j;
                if (des_pid == 0) {
                    copy_2d_array(block, tmp_block, bsz); 
                } 
                else {
                    // send block to other processes
                    MPI_Send(&tmp_block[0][0], bsz*bsz, MPI_DOUBLE, des_pid, des_pid, MPI_COMM_WORLD);
                }

                free_2d_double(tmp_block, bsz, bsz);
            }
        }        
        cout << "In process " << id << ", received block is: ";
        print_2d_array(block, bsz, bsz);
        free_2d_double(G, N, N);
    }
    else {
        MPI_Recv(&block[0][0], bsz*bsz, MPI_DOUBLE, 0, id, MPI_COMM_WORLD, &status);
        cout << "In process " << id << ", received block is: ";
        print_2d_array(block, bsz, bsz);
    }

    /***************************************************************************************************/
    // Spit into row- and col-wise groups
    int row_group_id = id / d;
    int col_group_id = (id % d) + d;

    MPI_Comm rowComm;
    MPI_Comm colComm;
    MPI_Comm_split(MPI_COMM_WORLD, row_group_id, id, &rowComm);
    MPI_Comm_split(MPI_COMM_WORLD, col_group_id, id, &colComm);
     
    int row_rank, col_rank;
    int row_group_sz, col_group_sz;
    MPI_Comm_size(rowComm, &row_group_sz);
    MPI_Comm_rank(rowComm, &row_rank);
    MPI_Comm_size(colComm, &col_group_sz);
    MPI_Comm_rank(colComm, &col_rank);

    MPI_Request row_group_forward, col_group_forward;
    MPI_Request left, right, up, down;
    // cout << "In process " << id << ", row_group/row_rank is: " << row_group_id << "/" << row_rank << endl;
    // cout << "In process " << id << ", col_group/col_rank is: " << col_group_id << "/" << col_rank << endl;

    double * kth_col = (double *) malloc(bsz * sizeof(double));
    double * kth_row = (double *) malloc(bsz * sizeof(double));

    for (int k = 0; k < N; k++) {
        // I have row data to send in kth iteration, 
        // if k are in the interval of [row_rank*bsz, (row_rank+1)*bsz]
        if (k >= row_rank*bsz && k < (row_rank+1)*bsz) {
            // Send kth row to neighbors
            for (int i = 0; i < bsz; i++) {
                kth_col[i] = block[i][k % bsz];
            }
            // cout << "In process id " << id << ", " << k <<"th_row array is: ";
            // print_1d_array(kth_row, bsz);

            if (row_rank > 0) {
                // send left
                MPI_Isend(&kth_col[0], bsz, MPI_DOUBLE, row_rank - 1, k, rowComm, &left);
                cout << "In process id " << id << ", " << "send col in iteration " << k << ":";
                print_1d_array(kth_col, bsz);
            }
            if (row_rank < d-1) {
                // send right
                MPI_Isend(&kth_col[0], bsz, MPI_DOUBLE, row_rank + 1, k, rowComm, &right);
                cout << "In process id " << id << ", " << "send col in iteration " << k << ":";
                print_1d_array(kth_col, bsz);
            }
        } 
        else {
            // Receive data from neighbours
            MPI_Recv(&kth_col[0], bsz, MPI_DOUBLE, MPI_ANY_SOURCE, k, rowComm, &status);
            cout << "In process id " << id << ", " << "received in col in iteration " << k << ":";
            print_1d_array(kth_col, bsz);

            if (status.MPI_SOURCE > row_rank && row_rank > 0) {
                // forward left
                MPI_Isend(&kth_col[0], bsz, MPI_DOUBLE, row_rank - 1, k, rowComm, &row_group_forward);
                cout << "In process id " << id << ", " << "forward col array in iteration " << k << ":";
                print_1d_array(kth_col, bsz);
            }
            else if (row_rank < d-1) {
                // forward right
                MPI_Isend(&kth_col[0], bsz, MPI_DOUBLE, row_rank + 1, k, rowComm, &row_group_forward);
                cout << "In process id " << id << ", " << "forward col array in iteration " << k << ":";
                print_1d_array(kth_col, bsz);
            }
        }

        // I have col data to send in kth iteration
        // if k are in the interval of [col_rank*bsz, (col_rank+1)*bsz]
        if (k >= col_rank*bsz && k < (col_rank+1)*bsz) {
            // Send kth col to neighbours
            for (int j = 0; j < bsz; j++) {
                kth_row[j] = block[k % bsz][j];
            }
            // cout << "In process id " << id << ", " << k <<"th_col array is: ";
            // print_1d_array(kth_col, bsz);

            if (col_rank > 0) {
                // send up
                MPI_Isend(&kth_row[0], bsz, MPI_DOUBLE, col_rank - 1, k, colComm, &up);
                cout << "In process id " << id << ", " << "send row array in iteration " << k << ":";
                print_1d_array(kth_row, bsz);
            }
            if (col_rank < d-1) {
                // send down
                MPI_Isend(&kth_row[0], bsz, MPI_DOUBLE, col_rank + 1, k, colComm, &down);
                cout << "In process id " << id << ", " << "send row array in iteration " << k << ":";
                print_1d_array(kth_row, bsz);
            }
        }
        else {
            // Receive data from neighbours
            MPI_Recv(&kth_row[0], bsz, MPI_DOUBLE, MPI_ANY_SOURCE, k, colComm, &status);
            cout << "In process id " << id << ", " << "received in row in iteration " << k << ":";
            print_1d_array(kth_row, bsz);

            if (status.MPI_SOURCE > col_rank && col_rank > 0) {
                // forward up
                MPI_Isend(&kth_row[0], bsz, MPI_DOUBLE, col_rank - 1, k, colComm, &col_group_forward);
                cout << "In process id " << id << ", " << "forward row array in iteration " << k << ":";
                print_1d_array(kth_row, bsz);
            }
            else if (col_rank < d-1) {
                // forward down
                MPI_Isend(&kth_row[0], bsz, MPI_DOUBLE, col_rank + 1, k, colComm, &col_group_forward);
                cout << "In process id " << id << ", " << "forward row array in iteration " << k << ":";
                print_1d_array(kth_row, bsz);
            }
        }

        // Compute
        for (int i = 0; i < bsz; i++) {
            for (int j = 0; j < bsz; j++) {
                double old_ik = kth_col[i];
                double old_kj = kth_row[j];
                block[i][j] = min(block[i][j], old_ik + old_kj);
            }
        }

        // Send Barrier
        if (k >= row_rank*bsz && k < (row_rank+1)*bsz) {
            if (row_rank > 0) {
                MPI_Wait(&left, MPI_STATUS_IGNORE);
            }
            if (row_rank < d - 1) {
                MPI_Wait(&right, MPI_STATUS_IGNORE);
            }
        }
        else {
            if (row_rank > 0 && row_rank < d - 1) {
                MPI_Wait(&row_group_forward, MPI_STATUS_IGNORE);
            }
        }

        if (k >= col_rank*bsz && k < (col_rank+1)*bsz) {
            if (col_rank > 0) {
                MPI_Wait(&up, MPI_STATUS_IGNORE);
            }
            if (col_rank < d - 1) {
                MPI_Wait(&down, MPI_STATUS_IGNORE);
            }
        }
        else {
            if (col_rank > 0 && col_rank < d - 1) {
                MPI_Wait(&col_group_forward, MPI_STATUS_IGNORE);
            }
        }
    }


    cout << "In process id " << id << ", final result array is: ";
    print_2d_array(block, bsz, bsz);

    free(kth_col);
    free(kth_row);

    MPI_Comm_free(&rowComm);
    MPI_Comm_free(&colComm);
    MPI_Finalize();
    return 0;
}



// allocate contiguous memory for a 2d array in size rows*cols
double** alloc_2d_double(int rows, int cols) {
    double * data = (double *) malloc(rows * cols * sizeof(double));
    double ** array = (double **) malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        array[i] = &(data[cols * i]);
    }
    return array;
}

void free_2d_double(double ** G, int rows, int cols) {
    free(G[0]);
    free(G);
}

// generate random graph represented by 2d array, if no edge, weight is infinity
void generate_graph(double ** G, int num_of_vertices) {
    double INF = numeric_limits<double>::infinity(); 
    // srand(time(0)); 

    for (int i = 0; i < num_of_vertices; i++) {
        G[i][i] = 0;
        for (int j = i+1; j < num_of_vertices; j++) {
            G[i][j] = rand() % 5;
            if ((G[i][j] - 0.0) < 1E-3) {
                G[i][j] = INF;
            }
        }
        for (int j = 0; j < i; j++) {
            G[i][j] = G[j][i];
        }
    } 
}

void copy_2d_array(double ** destination, double ** source, int num) {
    for (int i = 0; i < num; i++) {
        for (int j = 0; j < num; j++) {
            destination[i][j] = source[i][j];
        }
    }
}

// print 2d array
void print_2d_array(double ** D, int rows, int cols) {
    cout << endl;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            cout << setw(4) << D[i][j];
        }
        cout << endl;
    }
    cout << endl;
}

void print_1d_array(double * a, int len) {
    for (int i = 0; i < len; i++) {
        cout << setw(4) << a[i];
    }
    cout << endl;
}