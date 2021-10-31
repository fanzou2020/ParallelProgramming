# include <cstdlib>
# include <ctime>
# include <iomanip>
# include <iostream>
# include <mpi.h>


/*


 * +-----------+ +-------------------+ +-----------+
 * | Process 0 | |    Process 1      | | Process 2 |
 * +-+-------+-+ +-+-------+-------+-+ +-+-------+-+
 *   | Value |     | Value | Value |     | Value |
 *   |  100  |     |  101  |  102  |     |  103  |
 *   +-------+     +-------+-------+     +-------+
 *         |            |     |                |
 *         |            |     |                |
 *         |            |     |                |
 *         |            |     |                |
 *         |            |     |                |
 *         |            |     |                |
 *       +-----------------------------------------+
 *       |                Process 0                |
 *       +-----+-----+-----+-----+-----+-----+-----+
 *       | 100 |  0  | 101 | 102 |  0  |  0  | 103 |
 *       +-----+-----+-----+-----+-----+-----+-----+



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
			
			//Declare Receive Array
			int recvArr[7]={};
			
			 // Number of elements received from each process
			int counts[3]= {1,2,1}; 
			
			//Starting position of each chunk to be received
			int displacements[3]={6, 2,0};  
			
			//Declare Send Value
			int sendVal=100;
	
		/*MPI_Gatherv(
					const void *sendbuf, 
					int sendcount, 
					MPI_Datatype sendtype,
                	void *recvbuf, 
					const int *recvcounts, 
					const int *displs,
                	MPI_Datatype recvtype, 
					int root, 
					MPI_Comm comm)
			*/
			
			
			MPI_Gatherv(&sendVal,1,MPI_INT,recvArr,counts,displacements,MPI_INT,masterProcess,MPI_COMM_WORLD);
			
			
			///***MPI will wait for scatter on all processes to complete before proceeding to next line. So you can print/ do operations on your received array right after scatter.
			//Still Master Process 0 here.
			
			for(int i=0;i<7;i++){
				
				cout << "\nP" << id <<" :  "<<recvArr[i]<<endl;
				
			}
	
			break;
		}
	
	
		case 1: {
			
			//Declare send array
			
			int sendArr[2]={101,102};
			
			MPI_Gatherv(sendArr,2,MPI_INT,NULL,NULL,NULL,MPI_INT,masterProcess,MPI_COMM_WORLD);
		
			break;			
		}
			
			
		case 2:  {

			//Declare send value
			int sendVal=103;
			
			MPI_Gatherv(&sendVal,1,MPI_INT,NULL,NULL,NULL,MPI_INT,masterProcess,MPI_COMM_WORLD);
		
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

