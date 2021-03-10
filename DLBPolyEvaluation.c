/*
  Matthew Dupont
  COSC 6060 Exercise 2
  2021-03-12
*/

/*
   To do: Complete the program to output the combined answer for polynomial evaluation by
   Process 0 (Master). Compare the correctness with sequential code.
   
   Use chunk size (or K) terms in message passing instead of 1 term. Use K = 10, 20.
*/

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#include <mpi.h>
#include <getopt.h>

//#define MAX 50000
//#define MAX 10
#define COEFFICIENT 1
//#define K 10

// Use Pascal server for this assignment. Use Everest only when you can not access Pascal.
// mpicc DLBPolyEvaluation.c
// mpirun -np 3 ./a.out

// Two tags to distinguish work messages from process-termination messages
#define WORKTAG 1
#define DIETAG 2
#define FINALTAG 3
#define LENGTHTAG 4
#define INDIVIDUALTIMETAG 5

// MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
// MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)

enum verbosity { v_veryVerbose=20, v_verbose=10, v_descriptive=0, v_terse=-10, v_veryTerse=-20 };

double power(double x, int degree)
{     
  if(degree == 0)  return 1;
  
  if(degree == 1)  return x;

  return x * power(x, degree - 1);
}

// to be called from worker process
double evaluateTerm(int coefficient, int degree, double x, int verbosity)
{
  if (verbosity > v_verbose) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("Thread %.2i - degree %i coeff %i value %0.2f \n", rank, degree, coefficient, x);
    fflush(stdout);
  }
  return coefficient * power(x, degree);
}

double evaluateTerms(int* terms, int termLength, double x, int verbosity)
{
  int i;
  double localTotal = 0.0;
  for (i = 0; i < termLength; i++)
  {
    int degree = terms[2*i];
    int coefficient = terms[(2*i)+1];
    //printf("Thread: %i; Degree: %i; Coefficient:%i\n", rank, degree, coefficient);
    localTotal += evaluateTerm(coefficient, degree, x, verbosity);
  }

  if (verbosity > v_descriptive) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("Thread %.2i - Chunk of size %i produced total %f\n", rank, termLength, localTotal);
    fflush(stdout);
  }

  return localTotal;
}

double sequential(int coeffArr[], int numPolynomials, double x)
{
  //Sequential was not including the final term
  int i;
  double answer = 0;
  
  for( i = 0; i < numPolynomials;  i++)
  {
    double powerX = power(x, i);

    answer = answer + coeffArr[i] * powerX;
  }
  return answer;
}

void initialize(int coeffArr[], int numPolynomials)
{
   int i;
   for( i = 0; i < numPolynomials; i++)
   {
      coeffArr[i] = COEFFICIENT;
   }
}

void packAndSendTerm(int coeffArr[], int coeffLength, int* coeffPosition, int term[], int termLength, int targetWorker, int verbosity)
{
  if (verbosity > v_verbose) printf("Master - Packing and sending terms: position: %i terms: %i\n", *coeffPosition, termLength);

  int i;
  for (i = 0; i < termLength; i++){
    term[2*i] = *coeffPosition;               // degree
    term[2*i+1] = coeffArr[*coeffPosition];   // coefficient
    (*coeffPosition) ++;
    if (*coeffPosition >= coeffLength) break;
  }
  MPI_Send(term, 2 * termLength, MPI_INT, targetWorker, WORKTAG, MPI_COMM_WORLD); 
}

double queueMaster(int coeffArr[], int numPolynomials, int chunkSize, int verbosity)
{
  int numProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  MPI_Status status;
  int tCount = 0;
  int workerRank;
  int *term = (int *)malloc(sizeof(int) * 2 * chunkSize);
  double totalResult = 0.0;

  if (verbosity > v_terse) printf("Master - Seeding initial work...\n");
  //seed initial workloads
  for(workerRank = 1; workerRank < numProcs; workerRank++)
  { 
    packAndSendTerm(coeffArr, numPolynomials, &tCount, term, chunkSize, workerRank, verbosity);
  }

  if (verbosity > v_terse) printf("Master - Sending work from queue...\n");
  // receive the work done by workers
  double termAnswer;
  while(tCount < numPolynomials)
  {    
    // Receive the results from workers and accumulate the result
    MPI_Recv(&termAnswer, 1, MPI_DOUBLE, MPI_ANY_SOURCE, WORKTAG, MPI_COMM_WORLD, &status);
    totalResult += termAnswer;
    // Send new results to workers.
    workerRank = status.MPI_SOURCE;

    //If numPolynomials is not a multiple of chunkSize, we need to specially handle the last values.
    if (tCount > (numPolynomials - chunkSize)) {
      int remainingTerms = (numPolynomials - tCount);
      if (verbosity > v_terse) printf("Master - Position = %i; Sending %i final values...\n", tCount, remainingTerms);
      bzero(term, sizeof(int)*2*chunkSize);
      //Send dummy message with FINALTAG
      MPI_Send(term, 2* chunkSize, MPI_INT, workerRank, FINALTAG, MPI_COMM_WORLD);
      //Send length of final truncated message
      MPI_Send(&remainingTerms, 1, MPI_INT, workerRank, LENGTHTAG, MPI_COMM_WORLD);
      //Send final truncated term.
      packAndSendTerm(coeffArr, numPolynomials, &tCount, term, remainingTerms, workerRank, verbosity);
      if (verbosity > v_terse) printf("Master - Finished sending final values...\n");
    } else {
      packAndSendTerm(coeffArr, numPolynomials, &tCount, term, chunkSize, workerRank, verbosity);
    }
  }

  if (verbosity > v_terse) printf("Master - Sending die tags...\n");
  // Send (tag = DIETAG) to all workers
  for(workerRank = 1; workerRank < numProcs; workerRank++)
  {
    // Receive all outstanding work.
    MPI_Recv(&termAnswer, 1, MPI_DOUBLE, MPI_ANY_SOURCE, WORKTAG, MPI_COMM_WORLD, &status);
    totalResult += termAnswer;
    int worker = status.MPI_SOURCE;
    // Zero the term array - data in it is no longer needed.
    bzero(term, sizeof(int)*2*chunkSize);
    // Send die tag to worker
    MPI_Send(term, 2 * chunkSize, MPI_INT, worker, DIETAG, MPI_COMM_WORLD);
  }

  free(term);
  return totalResult;
}
  
void worker(double x, int chunkSize, int verbosity)
{

  // Start timer
  double elapsed_time = - MPI_Wtime();

	MPI_Status status;
  int val;
  int *term = (int *)malloc(sizeof(int) * 2 * chunkSize); 
        
  // worker keeps looping; waiting for work items from master process, 
  // unless it gets termination message
  while(1)
  { 
    MPI_Recv(term, 2*chunkSize, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    if(status.MPI_TAG == FINALTAG)
    {
      int remainingLength;
      MPI_Recv(&remainingLength, 1, MPI_INT, 0, LENGTHTAG, MPI_COMM_WORLD, &status);
      MPI_Recv(term, 2*remainingLength, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      double answer = evaluateTerms(term, remainingLength, x, verbosity);
      MPI_Send(&answer, 1, MPI_DOUBLE, 0, WORKTAG, MPI_COMM_WORLD);
    }
    if(status.MPI_TAG == DIETAG)
    {
      // End timer 
      elapsed_time = elapsed_time + MPI_Wtime();

      if (verbosity > v_descriptive) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        printf("Thread %.2i - TERMINATING AFTER %f seconds. BYE \n", rank, elapsed_time);
        fflush(stdout);
      }
      MPI_Send(&elapsed_time, 1, MPI_DOUBLE, 0, INDIVIDUALTIMETAG, MPI_COMM_WORLD);
      return;
    }
    else    // (status.MPI_TAG = WORKTAG)
    {
      double answer = evaluateTerms(term, chunkSize, x, verbosity);
      MPI_Send(&answer, 1, MPI_DOUBLE, 0, WORKTAG, MPI_COMM_WORLD);
    }
  } 
  free(term);
}
 

//Populates the provided buffer with timings from each proc, with index = procnum.
void aggregateDetailedTimings(int rank, int numProcs, double elapsed_time, double* buffer) {
  MPI_Status status;
  int i;
  if (rank == 0) {
    buffer[0] = elapsed_time;
    for (i = 1; i < numProcs; i++) {
      MPI_Recv(&buffer[i], 1, MPI_DOUBLE, i, INDIVIDUALTIMETAG, MPI_COMM_WORLD, &status);
    }
  }
}

// Driver Code
int main(int argc, char **argv)
{
  //Default values
  int chunkSize = 10;
  int numPolynomials = 50000;
  double variable = .99;

  int verbosity = v_descriptive;

  //Calculate Inputs
  static struct option longOpts[] =
  {
      // This does involve setting values in getOpt, which is thread unsafe, 
      // but MPI should be multi-process and thus this shouldn't be a concern.
      // See https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Options.html

      //long option,    argument status,    flag, short option
      { "help",	        no_argument,	      0,    'h' },
      { "poly_size",	  required_argument,	0,    'p' },
      { "variable",	    required_argument,  0,    'x' },
      { "chunksize",	  required_argument,  0,    'c' },
      { "veryVerbose",  no_argument,        0,    'w' },
      { "verbose",	    no_argument,        0,    'v' },
      { "descriptive",  no_argument,        0,    'd' },
      { "terse",	      no_argument,        0,    't' },
      { "veryTerse",    no_argument,        0,    'y' }
  };

  int c;
  int optionIndex = 0;
  int error = 0;

  while (error == 0 && (c = getopt_long(argc, argv, "hp:x:c:wvdty", longOpts, &optionIndex)) > -1)
  {
    switch(c)
    {
      //Usage
      case 'h': 
        printf("TODO\n");
        break;
      //Parameters
      case 'p':
        numPolynomials = atoi(optarg);
        break;
      case 'x':
        variable = atof(optarg);
        break;
      case 'c':
        chunkSize = atoi(optarg);
        break;
      //Logging
      case 'w':
        verbosity = v_veryVerbose;
        break;
      case 'v':
        verbosity = v_verbose;
        break;
      case 'd':
        verbosity = v_descriptive; 
        break;
      case 't':
        verbosity = v_terse; 
        break;
      case 'y':
        verbosity = v_veryTerse; 
        break;

      default:
        error = 1;
        break;
    }
  }  

  int *coeffArr = (int *)malloc(sizeof(int) * numPolynomials);
  double x = variable;
  
  initialize(coeffArr, numPolynomials);
  
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int numProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  if (verbosity > v_terse && rank == 0) printf("Evaluating a polynomial of %i terms, in chunks of %i elements, with %i processes\n", numPolynomials, chunkSize, numProcs);
  
  // Start timer
  double elapsed_time = - MPI_Wtime();

  double result;

  if (numProcs <= 1 && rank == 0) 
  {
    //Use SingleThread
    if (verbosity > v_terse) printf("Evaluating as sequential...\n");
    result = sequential(coeffArr, numPolynomials, x);
  }
  else 
  {
    //Use MultiThread
    if (rank == 0) 
    {
      if (verbosity > v_terse) printf("Evaluating with Queue...\n");
      result = queueMaster(coeffArr, numPolynomials, chunkSize, verbosity);
    } 
    else 
    {
      worker(x, chunkSize, verbosity);
    }
  }

  // End timer 
  elapsed_time = elapsed_time + MPI_Wtime();
  
  double* times = (double *)malloc(sizeof(double) * numProcs);
  bzero(times, sizeof(double) * numProcs);
  if(verbosity == v_veryTerse)
  {
    aggregateDetailedTimings(rank, numProcs, elapsed_time, times);
  }

  if (rank == 0)
  {
    switch (verbosity) {
      case v_terse: 
      {
        printf("%i,%i,%f,%f,%f\n", numProcs, numPolynomials, variable, result, elapsed_time);
        break;
      }
      case v_veryTerse: 
      {
        int i;
        for (i = 0; i < numProcs; i++) {
          printf("%i,%i,%f,%f,%i,%f\n", numProcs, numPolynomials, variable, result, i, times[i]);
        }
        break;
      }
      default:
      {
        printf("Total Result: %f Max Time: %f\n", result, elapsed_time);
        break;
      }
    }
  }

  free(coeffArr);
  free(times);
  
  MPI_Finalize();
    
	return 0;
}