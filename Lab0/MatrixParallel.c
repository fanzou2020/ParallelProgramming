#include <stdio.h>
#include <mpi.h>

#define N 1000                // matrix dimension N * 32
#define NUM_ELEMS N*32      // number of elements in matrix, N * 32 

/*
matrix A dimension is n x m, matrix B dimension m x p
return result dimension n x p
*/
void matrixMult(int* result, int* matrixA, int* matrixB, int n, int m, int p);

/*
print matrix into standard output
*/
void printMatrix(int* matrix, int row, int col);

/**
 * The strategy is dividing matrix A by row,
 * for example, if we have two processes,
 * divide A[10][32] -> A1[5][32] (process 1), A2[5][32] (process 0);
 * divide A[21][32] -> A1[10][32] (process 1), A2[11][32] (process 0).
 * Dispatch each sub-matrix of A to each process, 
 * and then merge them together.
 */
int main(int argc, char** argv) {
    int myRank;     // rank
    int NUM_PROCS;  // number of processes
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &NUM_PROCS);

    int NUM_ROWS_PER_PROC = N / NUM_PROCS;
    int NUM_ELEMS_PER_PROC = NUM_ROWS_PER_PROC * 32;

    int sendBuff[NUM_ELEMS];          // matrix A
    int recvBuff[NUM_ELEMS_PER_PROC]; // sub-part of matrix A

    // Construct matrix B, all process need the content of matrix B.
    int matrixB[NUM_ELEMS];
    for (int i = 0; i < NUM_ELEMS; i++) {
        matrixB[i] = i / 32;
    }


    if (myRank == 0) {
        // Master process
        double startTime = MPI_Wtime();

        // initialize matrix A to sendBuff
        for (int i = 0; i < NUM_ELEMS; i++) {
            sendBuff[i] = i / 32;
        }

        // Divide matrix A, send different sub-matrix to other processes
        int j = 0;
        for (j = 0; j < NUM_PROCS-1; j++) {
            int startIndex = j * NUM_ELEMS_PER_PROC;

            // send sub-matrix to other processes
            MPI_Send(sendBuff + startIndex, NUM_ELEMS_PER_PROC, MPI_INT, j+1, j+1, MPI_COMM_WORLD);
        }

        // final total result matrix
        int result[N * N];

        // calculate sub-matrix assigned to master process, 
        // if N is odd, the remaining extra row is assigned to this process
        int resultStartIndex = j * N * NUM_ROWS_PER_PROC;
        int matrixA_StartIndex = j * NUM_ELEMS_PER_PROC;
        int rowsLeft = N - j * NUM_ROWS_PER_PROC;
        matrixMult(result + resultStartIndex, sendBuff + matrixA_StartIndex, matrixB, rowsLeft, 32, N);

        /*
        printf("partial result from process 0, startIndex = %d\n", resultStartIndex);
        printMatrix(result + resultStartIndex, rowsLeft, N);
        printf("end of partial from process 0\n\n");
        */

        // put the particalResult returned by other process into the final total result
        for (int k = 0; k < NUM_PROCS-1; k++) {
            int beginIndex = k * N * NUM_ROWS_PER_PROC;
            MPI_Recv(result + beginIndex, N * NUM_ROWS_PER_PROC, MPI_INT, k+1, k+1, MPI_COMM_WORLD, &status);

            /*
            printf("Received result from process %d, startIndex = %d\n", k+1, beginIndex);
            printMatrix(result + beginIndex, NUM_ROWS_PER_PROC, N);
            printf("End of Received result\n\n");
            */

        }

        /*
        printMatrix(result, N, N);
        */

        double endTime = MPI_Wtime();

        printf("\nMatrix Multiplication (N = %d) took = %f on %d processes\n", N, (endTime - startTime), NUM_PROCS);

    }

    else {
        // Other process

        // receive submatrix of A
        MPI_Recv(&recvBuff, NUM_ELEMS_PER_PROC, MPI_INT, 0, myRank, MPI_COMM_WORLD, &status);

        int partialResult[NUM_ROWS_PER_PROC * N];

        // calculate partial result
        matrixMult(partialResult, recvBuff, matrixB, NUM_ROWS_PER_PROC, 32, N);

        /*
        printf("partial result from process %d\n", myRank);
        printMatrix(partialResult, NUM_ROWS_PER_PROC, N);
        printf("end of partial result from process %d\n\n", myRank);
        */

        // send partial result to master process
        MPI_Send(partialResult, N * NUM_ROWS_PER_PROC, MPI_INT, 0, myRank, MPI_COMM_WORLD);

    }

    MPI_Finalize();

    return 0;
}

/*
matrix A dimension is n x m, matrix B dimension m x p
return result dimension n x p
*/
void matrixMult(int* result, int* matrixA, int* matrixB, int n, int m, int p) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
            int sum = 0;
            for (int k = 0; k < m; k++) {
                // sum += A_ik * B_kj
                sum += matrixA[i*m+k] * matrixB[k*p+j];
            }
            // result_ij = sum
            result[i*p+j] = sum;
        }
    }
}

/*
print matrix into standard output
*/
void printMatrix(int* matrix, int row, int col) {
    for (int i = 0; i < row; i++) {
        printf("[ ");
        for (int j = 0; j < col-1; j++) {
            printf("%d, ", matrix[i*row + j]);
        }
        printf("%d ]\n", matrix[i*row + col-1]);
    }
}










