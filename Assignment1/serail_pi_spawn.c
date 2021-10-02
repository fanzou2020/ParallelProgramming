#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#define ROUND 100

double dboard(double threshold);

int main(int argc, char *argv[]) 
{
    /* Some declaration here */
    int n;   // the same computation amount as n works in parallel version

    if (argc != 2) {
        printf("Usage: ./<program name> n. The computational work is equal of n tasks of parallel version.\n");
    } else {
        n = atoi(argv[1]);
    }

    clock_t start_time = clock();

    int p;
    double pi_sum_over_proccess;
    for (p = 0; p < n; p++) {
        srandom(p);
        int i;
        double pi_sum = 0, pi;

        for (i = 0; i < ROUND; i++) {
            pi_sum += dboard(1E-5);
        }
        pi = pi_sum / ROUND;
        pi_sum_over_proccess += pi;
        printf("pi = %f\n", pi);
    }
    printf( "calculated pi value is %f\n", pi_sum_over_proccess/n );
    clock_t end_time = clock();
    double execution_time = (double) (end_time - start_time) / CLOCKS_PER_SEC;
    printf("Serial execution time is %f\n", execution_time);
    return 0;
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