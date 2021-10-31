# include <cstdlib>
# include <ctime>
# include <iomanip>
# include <iostream>
# include <mpi.h>

using namespace std;

int main ( int argc, char *argv[] ) {

  int id;  //process Id  (Rank in MPI)

  int numProcesses;

  double startTime, elapsedTime;

//  Initialize MPI.

  MPI_Init ( &argc, &argv );  

  startTime = MPI_Wtime ( );

//  Get the number of processes.

  MPI_Comm_size ( MPI_COMM_WORLD, &numProcesses );

//  Get the individual process ID.

  MPI_Comm_rank ( MPI_COMM_WORLD, &id );


 //Number of elements to send to each process (size of each chunk of array)
 int count=2; 

//Initial Send Array
 int sendArr[count*numProcesses];

// Initialise Receive Buffer on All processes
  int recvArr[numProcesses*count]={};
	
	//Populate Initial Send Array
	for(int i=0;i<numProcesses*count;i++){
			sendArr[i]=  (id*10) +  (i+1);
	}	
	
	
	/*
	
	PO: 1 2 3 4 5 6...
	P1: 11 12 13 14 15 16...
	P2: 21 22 23 24 25 26...
	
	*/

	/*MPI_Alltoall(void* buffer_send,
                 int count_send,
                 MPI_Datatype datatype_send,
                 void* buffer_recv,
                 int count_recv,
                 MPI_Datatype datatype_recv,
                 MPI_Comm communicator);*/

    MPI_Alltoall(sendArr, count, MPI_INT, recvArr, count, MPI_INT, MPI_COMM_WORLD);


//All processes will print the final result.

	for(int i=0;i<numProcesses*count;i++){
				
		cout << "P" << id <<" :  "<<recvArr[i]<<"    ";
	}
	cout<<endl;

//  Terminate MPI.

  elapsedTime = MPI_Wtime ( ) - startTime;

  MPI_Finalize ( );

//  Process 0 prints a termination message.

  if ( id == 0 )
  {
    cout << "\nP" << id << ":    Elapsed time = " << elapsedTime << " seconds.\n";	
    cout << "P" << id << ":    Normal end of execution.\n";
  }

  return 0;
}

