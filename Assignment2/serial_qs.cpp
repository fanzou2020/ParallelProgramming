#include <time.h>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"


// how many number to be sorted
#define N 10000

// max value of the numbers
#define MAX_VALUE 10000

using namespace std;


int compare (const void * a, const void * b) {
    return ( *(int*)a - *(int*)b );
}
void print_array(int a[], int len) {
    for (int i = 0; i < len; i++) {
        cout << a[i] << " ";
    }
    cout << endl;
}

int main() {
    double start_time = MPI_Wtime();
    int a[N];
    for (int i = 0; i < N; i++) {
        a[i] = rand() % MAX_VALUE; 
        // a[i] = i;
    }
    qsort(a, N, sizeof(int), compare);
    print_array(a, N);
    double end_time = MPI_Wtime();
    double execution_time = end_time - start_time;
    cout << "Parallel quick sort time is: " << execution_time << endl;


}

