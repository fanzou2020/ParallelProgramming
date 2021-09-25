#include <stdio.h>
#include <time.h>
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

int main() {
    int matrixA[NUM_ELEMS];
    int matrixB[NUM_ELEMS];

    clock_t startTime = clock();

    for (int i = 0; i < NUM_ELEMS; i++) {
        matrixA[i] = i / 32;
        matrixB[i] = i / 32;
    } 

    int result[N * N];

    matrixMult(result, matrixA, matrixB, N, 32, N);

    /*
    printMatrix(result, N, N);
    */

    clock_t endTime = clock();
    double time_spent = (double)(endTime - startTime) / CLOCKS_PER_SEC;

    printf("\nMatrix Multiplication (N = %d) Sequentially took = %f\n", N, time_spent);

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



