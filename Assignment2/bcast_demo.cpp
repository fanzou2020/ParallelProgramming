#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include "mpi.h"

using namespace std;


int main ( int argc, char* argv[]) {
    
    int id;

    int numProcesses;

    MPI_Init( &argc, &argv );

    MPI_Comm_size( MPI_COMM_WORLD, &numProcesses);

    MPI_Comm_rank( MPI_COMM_WORLD, &id );

    int data = id * 10;

    MPI_Bcast( &data, 1, MPI_INT, 0, MPI_COMM_WORLD );

    cout << "in process " << id << " data is " << data << endl;

    MPI_Finalize(); 

}

