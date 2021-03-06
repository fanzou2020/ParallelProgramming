#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#define ROUND 100

// Keeps calculating pi until it converges
double dboard(double threshold);

int main(int argc, char *argv[])
{
    int processId, numChildren;

    double pi, sum;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &numChildren);

    MPI_Comm_rank(MPI_COMM_WORLD, &processId);

    MPI_Comm commParent;

    MPI_Comm_get_parent(&commParent);

    // calculate pi 
    srandom(processId);
    int i;
    double pi_sum;

    for (i = 0; i < ROUND; i++) {
        pi_sum += dboard(1E-5);
    }
    pi = pi_sum / ROUND;
    printf("calculate pi in %d, pi = %f\n", processId, pi);

    MPI_Reduce(&pi, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, commParent);

    MPI_Finalize(); 
    
}

double dboard(double threshold) {
    #define sqr(x)  ((x) * (x))
    long random(void);
    double x_coord, y_coord, pi, r;
    int score, n;
    unsigned int cconst;  /* must be 4-bytes in size */

    if (sizeof(cconst) != 4) {
        printf("Wrong data size for cconst variable in dboard routine!\n");
        printf("See comments in source file. Quitting.\n");
        exit(1);
    }

    cconst = 2 << (31 - 1);
    score = 0;

    /*throw darts at board until converges */
    double pi_new;
    n = 0;
    while (1) {
        r = (double) random()/cconst;
        x_coord = (2.0 * r) - 1.0;
        r = (double) random()/cconst;
        y_coord = (2.0 * r) - 1.0;
        
        // if dart lands in circle, increment score
        if ((sqr(x_coord) + sqr(y_coord)) <= 1.0) {
            score++;
        }
        pi_new = 4.0 * (double)score/(double)(++n);
        // printf("%d step, pi_new = %lf, pi = %lf\n", n, pi, pi_new);

        if ((fabs(pi_new - pi) < threshold) && n > 10000) 
            break; 
        else 
            pi = pi_new;
    }
    return pi;
}