# include <cstdlib>
# include <ctime>
# include <iomanip>
# include <iostream>
# include <mpi.h>


/*

 *
 *       +-----------------------------------------+
 *       |                Process 0                |
 *       +-----+-----+-----+-----+-----+-----+-----+
 *       | 100 |  0  | 101 | 102 |  0  |  0  | 103 |
 *       +-----+-----+-----+-----+-----+-----+-----+
 *         |            |     |                |
 *         |            |     |                |
 *         |            |     |                |
 *         |            |     |                |
 *         |            |     |                |
 *         |            |     |                |
 * +-----------+ +-------------------+ +-----------+
 * | Process 0 | |    Process 1      | | Process 2 |
 * +-+-------+-+ +-+-------+-------+-+ +-+-------+-+
 *   | Value |     | Value | Value |     | Value |
 *   |  100  |     |  101  |  102  |     |  103  |
 *   +-------+     +-------+-------+     +-------+


*/


using namespace std;

int main ( int argc, char *argv[] ) {

  int id;  //process Id  (Rank in MPI)

  int numProcesses;

  double startTime, elapsedTime;

  int masterProcess=0; //P0 will populate main array and scatter

//  Initialize MPI.

  MPI_Init ( &argc, &argv );  

  startTime = MPI_Wtime ( );

//  Get the number of processes.

  MPI_Comm_size ( MPI_COMM_WORLD, &numProcesses );

  if(numProcesses!=3){
  
     cout<<"This application is meant to be run with 3 processes.\n";
	 
     MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  

//  Get the individual process ID.

  MPI_Comm_rank ( MPI_COMM_WORLD, &id );


//Process 0 populates array and scatters it to others

switch(id) {
		
		case 0:  {
			
			//Populate Main Array
			int mainArr[7]={100,0,101,102,0,0,103};  
			
			 // Number of elements I need to send to each process
			int counts[3]= {1,2,1}; 
			
			//Starting position of each chunk or array to be sent
			int displacements[3]={6, 2,0};  
			
			//Declare Receive Value
			int recvVal;
	
	        			/* MPI_Scatterv(const void* buffer_send,
                 const int counts_send[],
                 const int displacements[],
                 MPI_Datatype datatype_send,
                 void* buffer_recv,
                 int count_recv,
                 MPI_Datatype datatype_recv,
                 int root,
                 MPI_Comm communicator);
			*/
			
			
			MPI_Scatterv(mainArr,counts,displacements,MPI_INT,&recvVal,1,MPI_INT,masterProcess,MPI_COMM_WORLD);
			
			cout << "\nP" << id <<" :  "<<recvVal<<endl;
	
			break;
		}
	
	
		case 1: {
			
			//Declare receive buffer
			int recvArr[2];
			
			MPI_Scatterv(NULL,NULL,NULL,MPI_INT,recvArr,2,MPI_INT,masterProcess,MPI_COMM_WORLD);
			
			cout << "\nP" << id <<" :  "<<recvArr[0]<<",   "<<recvArr[1]<<endl;
		
			break;			
		}
			
			
		case 2:  {

			//Declare receive value
			int recvVal;
			
			MPI_Scatterv(NULL,NULL,NULL,MPI_INT,&recvVal,1,MPI_INT,masterProcess,MPI_COMM_WORLD);
			
			cout << "\nP" << id <<" :  "<<recvVal<<endl;
		
			break;
		}
		
}  //End Of Switch


//  Terminate MPI.

  elapsedTime = MPI_Wtime ( ) - startTime;

  MPI_Finalize ( );

//  Process 0 prints a termination message.

  if ( id == 0 )
  {
    cout << "\n\nP" << id << ":    Elapsed time = " << elapsedTime << " seconds.\n";	
    cout << "P" << id << ":    Normal end of execution.\n";
  }

  return 0;
}

