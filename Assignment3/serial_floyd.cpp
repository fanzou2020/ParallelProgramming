#include <iostream>
#include <algorithm>
#include <limits>
#include <stdlib.h>

using namespace std;


// 2D array of size N*N
double** floyd_all_pair(int N, double ** G) {
    double D_old[N][N];  // allocate D_old in stack, initialize with G
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            D_old[i][j] = G[i][j];
        }
    }
    // allocate a new D[N][N] array in heap
    double** D_new = (double **) calloc(N, sizeof(double*));
    for (int i = 0; i < N; i++) {
        D_new[i] = (double *) calloc(N, sizeof(double));
    }

    for (int k = 0; k < N; ++k) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                D_new[i][j] = min(D_old[i][j], D_old[i][k] + D_old[k][j]);
            }
        }

        // copy D_new to D_old
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                D_old[i][j] = D_new[i][j];
            }
        }
    }
    return D_new;
}

void print_2d_array(double ** D, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << D[i][j] << " " ;
        }
        cout << endl;
    }
}



int main() {

    double INF = numeric_limits<double>::infinity();
    int N = 4;
    double G_2d[4][4] = {{0, 2, 3, 1}, {2, 0, INF, INF}, {3, INF, 0, 3}, {1, INF, 3, 0}};
    double * G[N];
    for (int i = 0; i < N; i++) {
        G[i] = G_2d[i];
    }

    double** res = floyd_all_pair(4, G);

    print_2d_array(res, N);

    return 0;

}
















