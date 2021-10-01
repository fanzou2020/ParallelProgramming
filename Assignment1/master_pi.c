#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

#define ROUND 2

int main(int argc, char *argv[])
{
    int numChildren, processId, numParents;

    double pi_sum;

    MPI_Comm childComm;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &numParents);

    MPI_Comm_rank(MPI_COMM_WORLD, &processId);

    printf("Parent task id is: %d. Total tasks in parent's world are: %d \n", processId, numParents);    
    double start_time = MPI_Wtime();

    if (argc != 2) {
        printf("usage: %s <number of worker>\n", argv[0]);
    }
    else {
        numChildren = atoi( argv[1] );  // Given by user in command line, 'mpirun -n 1 master 4'
    }

    MPI_Comm_spawn("child_pi", argv, numChildren, MPI_INFO_NULL, 0, MPI_COMM_SELF, &childComm, MPI_ERRCODES_IGNORE); 

    MPI_Reduce(&processId, &pi_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_ROOT, childComm);

    double pi_avg = pi_sum / numChildren;
    printf("pi_avg = %f, number of children is %d\n", pi_avg, numChildren);

    double end_time = MPI_Wtime();
    double execution_time = end_time - start_time;
    printf("Execution time using MPI_Spawn is %f\n", execution_time);

    MPI_Finalize();
    
}
