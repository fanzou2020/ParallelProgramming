# include <cstdlib>
# include <ctime>
# include <iomanip>
# include <iostream>
# include "mpi.h"

using namespace std;

int main ( int argc, char *argv[] ) {

  int id;  //process Id  (Rank in MPI)

  int numProcesses;

  //double startTime, elapsedTime;

//  Initialize MPI.

  MPI_Init ( &argc, &argv );  

  //startTime = MPI_Wtime ( );

//  Get the number of processes.

  MPI_Comm_size ( MPI_COMM_WORLD, &numProcesses );

//  Get the individual process ID.

  MPI_Comm_rank ( MPI_COMM_WORLD, &id );
  
  //I want 4 rows here...You can change
  int numRows=4;

//Group processes into rows
int color=id/numRows;  

//New communicator for each row
MPI_Comm rowComm;

/*MPI_Comm_split(MPI_Comm old_communicator,
                   int colour,
                   int key,
                   MPI_Comm* new_communicator);
	*/

//Color determines to which group my process goes to. Processes with same color go into the same group
//Key determines the rank of the proces in the new group. Lower the key, Lower the rank.
//** If 2 processes in same group have same key, the process with the initial lower rank gets the lower rank in the new group

MPI_Comm_split(MPI_COMM_WORLD, color, id, &rowComm);


//ID OF PROCESS IN NEW GROUP 
int rowId;

// SIZE OF NEW GROUP
int rowSize;  

MPI_Comm_size(rowComm, &rowSize);
MPI_Comm_rank(rowComm, &rowId);


cout<<"\nWorldRank/Size: "<<id<<"/"<<numProcesses<<"   RowRank/Size: "<<rowId<<"/"<<rowSize<<"  in Group: "<<color<<endl;


//Free your communicators after use
MPI_Comm_free(&rowComm);


//  Terminate MPI.

  //elapsedTime = MPI_Wtime ( ) - startTime;

  MPI_Finalize ( );

  return 0;
}

