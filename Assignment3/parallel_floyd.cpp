#include <iostream>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <algorithm>
#include "mpi.h"

// number of vertices
#define N 36

using namespace std;

double** alloc_2d_double(int rows, int cols);

void free_2d_double(double ** G, int rows, int cols);

void generate_graph(double ** G, int num_of_vertices);

void copy_2d_array(double ** destination, double ** source, int num);

void print_2d_array(double ** D, int rows, int cols);

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


    //***********************************************************************************
    // Split into row- and col-wise groups
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
    cout << "In process " << id << ", row_group/row_rank is: " << row_group_id << "/" << row_rank << endl;
    cout << "In process " << id << ", col_group/col_rank is: " << col_group_id << "/" << col_rank << endl;


    // Allocate 2d array to receive data from other processes in the group
    // [blocck blcok ... block] -> [block; block; ... block] in recved_rows
    double ** recved_rows = alloc_2d_double(N, bsz);
    double ** recved_cols = alloc_2d_double(N, bsz);

    // Broadcast data to other process in row-wise
    MPI_Bcast(&block[0][0], bsz*bsz, MPI_DOUBLE, row_rank, rowComm);
    // Receive data from other processes
    for (int i = 0; i < row_group_sz; i++) {
        if (i != row_rank) {
            MPI_Bcast(&recved_rows[i*bsz][0], bsz*bsz, MPI_DOUBLE, i, rowComm);
        } 
    }
    memcpy(&recved_rows[row_rank*bsz][0], &block[0][0], bsz*bsz*sizeof(double));

    /*
    if (row_rank == 0) {
        cout << "In process/row_rank " << id << "/" << row_rank << ", recved_rows array is:";
        print_2d_array(recved_rows, N, bsz);
    }
    */

    
    // Broadcast data to other process in col-wise
    MPI_Bcast(&block[0][0], bsz*bsz, MPI_DOUBLE, col_rank, colComm);
    // Receive data from other processes
    for (int i = 0; i < col_group_sz; i++) {
        if (i != col_rank) {
            MPI_Bcast(&recved_cols[i*bsz][0], bsz*bsz, MPI_DOUBLE, i, colComm);
        }
    }
    memcpy(&recved_cols[col_rank*bsz][0], &block[0][0], bsz*bsz*sizeof(double));
    
    /*
    if (col_rank == 0) {
        cout << "In process/col_rank " << id << "/" << col_rank << ", recved_cols array is:";
        print_2d_array(recved_cols, N, bsz);
    }
    */

    /*************************************************************************************/
    for (int k = 0; k < N; k++) {
        for (int i = 0; i < bsz; i++) {
            for (int j = 0; j < bsz; j++) {
                double old_ik = recved_rows[i + (k/bsz)*bsz][k % bsz];
                double old_kj = recved_cols[k][j];
                double tmp = min(block[i][j], old_ik + old_kj);
                block[i][j] = tmp;
            }
        }

        // Broadcast new block to other processes and construct new recved_rows, recved_cols
        MPI_Bcast(&block[0][0], bsz*bsz, MPI_DOUBLE, row_rank, rowComm);
        for (int r = 0; r < row_group_sz; r++) {
            if (r != row_rank) {
                MPI_Bcast(&recved_rows[r*bsz][0], bsz*bsz, MPI_DOUBLE, r, rowComm);
            }
        } 
        memcpy(&recved_rows[row_rank*bsz][0], &block[0][0], bsz*bsz*sizeof(double));

        MPI_Bcast(&block[0][0], bsz*bsz, MPI_DOUBLE, col_rank, colComm);
        for (int c = 0; c < col_group_sz; c++) {
            if (c != col_rank) {
                MPI_Bcast(&recved_cols[c*bsz][0], bsz*bsz, MPI_DOUBLE, c, colComm);
            }
        }
        memcpy(&recved_cols[col_rank*bsz][0], &block[0][0], bsz*bsz*sizeof(double));
    }

    cout << "In process id " << id << ", final result array is: ";
    print_2d_array(block, bsz, bsz);

    /**************************************************************************************/
    // Gather all blocks into one whole matrix result
    if (id == 0) {
        double ** shortest_path_matrix = alloc_2d_double(N, N);

        // root process receive result from other processes
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                double ** tmp_block = alloc_2d_double(bsz, bsz);

                int source_pid = i * d + j;
                if (source_pid == 0) {
                    copy_2d_array(tmp_block, block, bsz); 
                }
                else {
                    MPI_Recv(&tmp_block[0][0], bsz*bsz, MPI_DOUBLE, source_pid, source_pid, MPI_COMM_WORLD, &status);
                }

                // merge blocks into a whole graph matrix
                for (int k = 0; k < bsz; k++) {
                    for (int l = 0; l < bsz; l++) {
                        shortest_path_matrix[k + i * bsz][l + j * bsz] = tmp_block[k][l];
                    }
                }
                free_2d_double(tmp_block, bsz, bsz);
            }
        }
        cout << "The final shortest path matrix is :";
        print_2d_array(shortest_path_matrix, N, N);
        free_2d_double(shortest_path_matrix, N, N);

    } 
    else {
        MPI_Send(&block[0][0], bsz*bsz, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
    }


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

    ifstream file("demo_data/D0.txt");
    string str;
    int i = 0;
    while (getline(file, str)) {
        char * str_char_array = &str[0];
        // cout << str_char_array << endl;
        char * token = strtok(str_char_array, "  ");
        int j = 0;
        while (token != NULL) {
            // double d = stod(token);
            // printf("%s", token);
            int d = atoi(token);
            if (d == 0) G[i][j] = INF;
            else G[i][j] = double(d);
            if (i == j) G[i][j] = 0;
            token = strtok(NULL, "  ");
            j++;
        }
        i++;
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