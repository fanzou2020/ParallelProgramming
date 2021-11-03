#include <time.h>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

// how many number to be sorted
#define N 100

// max value of the numbers
#define MAX_VALUE 1000

using namespace std;

bool is_power_of_2(int x);

int partition(int a[], int len, int pivot);

void print_array(int a[], int len);

int compare (const void * a, const void * b) {
    return ( *(int*)a - *(int*)b );
}


int main( int argc, char* argv[] ) {

    int id;        // process Id (Rank in MPI) 

    int numProcesses;

    int masterProcess = 0;

    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses); 

    // Get the individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // if number of processes is not a power of 2, then exit program and print message
    // we want to simulate the hypercube topology in cluster
    if (!is_power_of_2(numProcesses)) {
        cout << "Pleasel use power of 2 process numbers" << endl;
        return 1;
    }

    int d = log2(numProcesses);    
    // cout << d << endl;

    /**************************** Scatter array to processes *************************/
    // calculate how many numbers in each process
    int count;
    if (N % numProcesses == 0) {
        count = N / numProcesses; 
    } else {
        count = N / numProcesses + 1;
    }

    // construct counts and displacements array for ScatterV 
    int counts[numProcesses];
    int displacements[numProcesses];

    for (int i = 0; i < numProcesses; i++) {
        displacements[i] = i * count;
        if (i != numProcesses - 1) {
            counts[i] = count;
        } else {
            counts[i] = N - count * (numProcesses-1);
        }
    }

    // receive buffer for each process
    int recv[counts[id]];

    // root process generate N numbers and scatter them to numProcesses processes
    if (id == 0) {

        // for (int j = 0; j < numProcesses; j++) {
        //     cout << "counts " << j << " = " << counts[j] << endl;
        //     cout << "displacements " << j << " = " << displacements[j] << endl;
        // }

        // genereate array with N random numbers
        int a[N];
        for (int i = 0; i < N; i++) {
            a[i] = rand() % MAX_VALUE; 
            // a[i] = i;
        }

        MPI_Scatterv(a, counts, displacements, MPI_INT, &recv, counts[id], MPI_INT, masterProcess, MPI_COMM_WORLD);

        cout << "process " << id << " received : "; 
        for (int k = 0; k < counts[id]; k++) {
            cout << recv[k] << " ";
        }
        cout << endl;
        
    }
    // other processes
    else {
        MPI_Scatterv(NULL, NULL, NULL, MPI_INT, &recv, counts[id], MPI_INT, masterProcess, MPI_COMM_WORLD);
        cout << "process " << id << " received : "; 
        for (int k = 0; k < counts[id]; k++) {
            cout << recv[k] << " ";
        }
        cout << endl;
    }
    


    /************************** Parallel Quicksort algorithm ************************/
    // Initialize buffer with data in recv[]
    int * buffer;   
    int len_buffer = counts[id];     
    buffer = (int *) malloc(len_buffer * sizeof(int));
    for (int i = 0; i < len_buffer; i++) {
        buffer[i] = recv[i];
    }
    
    MPI_Status status;

    for (int i = d; i > 0; i--) {
        // cout << "counts array is ";
        // print_array(counts, numProcesses);

        int color = id / pow(2, i);  // group id
        
        MPI_Comm rowComm;

        MPI_Comm_split(MPI_COMM_WORLD, color, id, &rowComm);

        // Id and size of process in new group
        int rowId, rowSize;
        MPI_Comm_size(rowComm, &rowSize);
        MPI_Comm_rank(rowComm, &rowId);

        cout << "\nWorldRank/size is " << id <<"/" << numProcesses << ", Rowrank/size is " << rowId << "/" << rowSize << " in Group:"
        << color << endl; 


        // choose pivot, and the root process in group broadcast pivot
        int pivot;

        if (rowId == 0) {
            int groupNum = numProcesses / pow(2, i);

            pivot = (MAX_VALUE / groupNum) * color + (MAX_VALUE / (groupNum * 2)) ; 
            cout << "pivot chosen is " << pivot << " in group " << color << endl;

            // Broadcast pivot in group
            MPI_Bcast(&pivot, 1, MPI_INT, 0, rowComm);

        } else {
            // Receive pivot from broadcast
            MPI_Bcast(&pivot, 1, MPI_INT, 0, rowComm);
        }


        // partition array buffer[] into B1 and B2, such that B1 <= pivot < B2
        cout << "buffer content is: ";
        print_array(buffer, len_buffer);
        cout << "len of buffer is " << len_buffer << endl;
        int startIndexB2 = partition(buffer, len_buffer, pivot);
        int* B1 = buffer;
        int* B2 = buffer + startIndexB2;
        int len_B1 = startIndexB2;
        int len_B2 = len_buffer - startIndexB2; 
        cout << "lenB1: " << len_B1 << ", and len_B2 " << len_B2 << endl;
        
        cout << "Process " << id << " after partition is: ";
        print_array(buffer, len_buffer);
        cout << "In rowId " << rowId << ", start index of B2 is " << startIndexB2 << endl;


        // find pair in group
        int pairId;
        if (rowId >= rowSize / 2) {
            pairId = rowId - rowSize / 2; 
            // cout << "rowId: " << rowId << " pairId:" << pairId << endl;

            // send B1 to pairId in group
            MPI_Send(&len_B1, 1, MPI_INT, pairId, 0, rowComm);  // send length of B1 with tag 0
            MPI_Send(B1, len_B1, MPI_INT, pairId, 1, rowComm);  // send B1 with tag 1

            // recv C from pairId
            int len_C;
            MPI_Recv(&len_C, 1, MPI_INT, pairId, 0, rowComm, &status);
            int C[len_C];
            MPI_Recv(C, len_C, MPI_INT, pairId, 1, rowComm, &status);
            cout << "in rowId " << rowId << ", C = ";
            print_array(C, len_C);

            // union B2 and C
            int len_new_buffer = len_B2 + len_C;
            int * new_buffer = (int *) malloc(len_new_buffer * sizeof(int));

            for (int i = 0; i < len_new_buffer; i++) {
                if (i < len_B2) {
                    new_buffer[i] = B2[i];
                } else {
                    new_buffer[i] = C[i - len_B2];
                }
            } 
            free(buffer);
            buffer = new_buffer;
            len_buffer = len_new_buffer;
            counts[id] = len_buffer;
            cout << "new buffer in pId/rowId " << id << "/" << rowId << " is:";
            print_array(buffer, len_buffer);

        
        } else {
            pairId = rowId + rowSize / 2;
            // cout << "rowId: " << rowId << " pairId:" << pairId << endl;

            // send B2 to pairId in group
            MPI_Send(&len_B2, 1, MPI_INT, pairId, 0, rowComm);  // send length of B2 with tag 0
            MPI_Send(B2, len_B2, MPI_INT, pairId, 1, rowComm);  // send B2 with tag 1

            // recv C from pairId 
            int len_C;
            MPI_Recv(&len_C, 1, MPI_INT, pairId, 0, rowComm, &status);
            int C[len_C];
            MPI_Recv(C, len_C, MPI_INT, pairId, 1, rowComm, &status);
            // cout << "in rowId " << rowId << ", C = " << " ";
            // print_array(C, len_C);

            // union B1 and C
            int len_new_buffer = len_B1 + len_C;
            int * new_buffer = (int *) malloc (len_new_buffer * sizeof(int));

            for (int i = 0; i < len_new_buffer; i++) {
                if (i < len_B1) {
                    new_buffer[i] = B1[i];
                } else {
                    new_buffer[i] = C[i - len_B1];
                }
            } 
            free(buffer);
            buffer = new_buffer;
            len_buffer = len_new_buffer;
            counts[id] = len_buffer;
            cout << "new buffer in pId/rowId " << id << "/" << rowId << " is:";
            print_array(buffer, len_buffer);

        }

        MPI_Comm_free(&rowComm);
    }


    /**************** local quicksort in each process **************/
    qsort(buffer, len_buffer, sizeof(int), compare);
    cout << "Final sorted buffer in pId " << id << " is:";
    print_array(buffer, len_buffer);


    /**************** Gather all elements in the order of process id, output the final sorted array **************/ 
    // each process other than root will send their buffer length to root process, and send buffer content by GatherV
    if (id == 0) {
        // construct counts and displacements array for GatherV
        int gather_counts[numProcesses];
        int gather_displacements[numProcesses];

        gather_counts[0] = len_buffer;
        for (int p = 1; p < numProcesses; p++) {
            MPI_Recv(&gather_counts[p], 1, MPI_INT, p, p, MPI_COMM_WORLD, &status);
        }
        cout << "gather counts array is ";
        print_array(gather_counts, numProcesses);

        gather_displacements[0] = 0;
        for (int d = 1; d < numProcesses; d++) {
            gather_displacements[d] = gather_displacements[d-1] + gather_counts[d-1];
        }
        cout << "gather displacements array is ";
        print_array(gather_displacements, numProcesses);

        // Receive buffern content
        int sorted_recv_array[N];
        MPI_Gatherv(buffer, len_buffer, MPI_INT, sorted_recv_array, gather_counts, 
                    gather_displacements, MPI_INT, 0, MPI_COMM_WORLD);

        cout << "\n\nThe Final Sorted Array is ";
        print_array(sorted_recv_array, N);

    } 
    else {
        // send buffer length
        MPI_Send(&len_buffer, 1, MPI_INT, 0, id, MPI_COMM_WORLD);

        // send buffer content
        MPI_Gatherv(buffer, len_buffer, MPI_INT, NULL, NULL, NULL, MPI_INT, 0, MPI_COMM_WORLD);

    }

    MPI_Finalize();
    
    return 0;

}


bool is_power_of_2(int x) {
    return x > 0 && !(x & (x-1));
}


// partition a[len] into a[0, 1, ... j-1] <= pivot < a[j, ... len-1], return j
int partition(int a[], int len, int pivot) {
    if (len == 0) {
        return 0;
    }

    // cout << "In partition function, pivot = " << pivot << endl;
    int low = 0, hi = len-1;
    
    while(1) {
        while (a[low] <= pivot) {
            low++;
            if (low == len) break;
        }
        while (a[hi] > pivot) {
            hi--;
            if (hi < 0) break;
        }

        if (low >= hi) break;
        
        int tmp = a[low];
        a[low] = a[hi];
        a[hi] = tmp;
    }
    return low;
}

void print_array(int a[], int len) {
    for (int i = 0; i < len; i++) {
        cout << a[i] << " ";
    }
    cout << endl;
}

