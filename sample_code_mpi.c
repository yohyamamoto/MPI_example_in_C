/*
	This program uses a random number generator and the Box-Muller (polar form) 
	transformation to generate normally distributed random numbers.
	I used mpi to demonstrate some for loops can get a benefit of parallel 
	computings (although this is probably not the best place to implement it).

	Compile this with: mpicc sample_code_mpi.c
	and run it with: mpirun ./a.out -np 4 
	or: mpirun -hostfile machinefile -np 4 ./a.out

	© March 19th, 2014	Yoh Yamamoto
	
*/
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/* function declaration */
void BoxMullar(double m, double s, double* output1, double* output2);
void DrawHistgram(double a[], int N);

int main(int argc, char *argv[]) {
	int i,size;
	int ierr, num_procs, my_id;	// MPI 
	int N=1000000;		// Number of random numbers to generate
	double a[N];
	double *a_global=NULL;	// An array to store results
	//double ran_num
	double output1, output2;	//Outputs from the Box-Muller transformation function
	double mean = 0.5;		// Mean of the Normal distribution
	double stddev = 0.1;	// Standard deviation of the Normall distribution
	time_t t;

	// Initialize MPI
	ierr = MPI_Init(&argc, &argv);

	// Get MPI rank and number of procs
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
        ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	/* Intializes random number generator. Here, I'm using time*rank as a random seed
	for a simpilicity sake. */
	srand((unsigned) time(&t)*my_id);	

	/* Before the Box-Muller transformation */
	size =  num_procs * N;
	if( my_id == 0 ) {
		a_global = malloc(size * sizeof(double) );
		printf("This is uniformly distributed random numbers from rand()\n");
	}
	for(i=0; i<N;++i){
		// Generate a number between 0 and 1
		a[i]=  (double) rand() / (RAND_MAX);
	}
	MPI_Gather( &a, N, MPI_DOUBLE, a_global, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if( my_id == 0 ) {
		DrawHistgram(a_global,size);
	}

	printf("hello from %d\n", my_id);
	free(a_global);

	/* After the Box-Muller transformation */

        if( my_id == 0 ) {
                a_global = malloc(size * sizeof(double) );
		printf( "This is normally distributed random numbers (Gaussian)\n");
	}
	for(i=0; i< (int) (N/2); ++i){
		BoxMullar(mean, stddev, &output1, &output2);
		//printf( "%20.18f\n%20.18f\n", output1, output2);
		a[2*i] = output1;
		a[2*i+1] = output2;
	}
	MPI_Gather( &a, N, MPI_DOUBLE, a_global, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if( my_id == 0 ) {
		DrawHistgram(a_global,num_procs * N);
	}

	ierr = MPI_Finalize();
	return 0;
}

/* The Box-Muller transformation 
	This function transforms the uniformly generated random numbers from rand()
	into normally distributed random numbers using Box-Muller transformations.

	This function takes mean and standard deviation and output two numbers.

	Math behind the transformation: 
	http://en.wikipedia.org/wiki/Box-Muller_transform
 */
void BoxMullar(double m, double s, double* output1, double* output2) {
	/* local variables */
	double x1,x2, y1, y2;
	double temp1, temp2;	// temporary  variables
	double w=2.0; // initial value
	//time_t t;
	//srand((unsigned) time(&t)*time(&t));
	while(1) {
		while(w >= 1.0) {
			x1=(2.0* (double) rand() / (RAND_MAX))-1.0;
			x2=(2.0* (double) rand() / (RAND_MAX))-1.0;
			w=x1*x1+x2*x2;
		}
		w=sqrt((-2.0*log(w))/w);
		y1=x1*w;
		y2=x2*w;
		temp1=m+y1*s;
		temp2=m+y2*s;
		if(temp1 > 0.0 && temp1 < 2.0*m) {
			if(temp2 > 0.0 && temp2 < 2.0*m) {
				break;
			}
		}
	}

	// Outputs: transformed random numbers.
	*output1 = temp1;
	*output2 = temp2;
}

/* A simple function to draw a histgram from a data set.
I would normally use gnuplot, matplotlib, or paraview for creating 
publication quality plots and figures.

   a[] is the data set and N is the array size
*/
void DrawHistgram(double a[], int N){
	int i,j;
	int hist[10];
	for(i=0;i<10;i++){
		hist[i] = 0;
	}
	// histgram bins
	for(i=0; i<N; i++) {
		if(a[i] < 0.1){
			hist[0]++;
		}
		else if (a[i] < 0.2) {
			hist[1]++;
		}
		else if (a[i] < 0.3) {
			hist[2]++;
		}
		else if (a[i] < 0.4) {
			hist[3]++;
		}
		else if (a[i] < 0.5) {
			hist[4]++;
		}
		else if (a[i] < 0.6) {
			hist[5]++;
		}
		else if (a[i] < 0.7) {
			hist[6]++;
		}
		else if (a[i] < 0.8) {
			hist[7]++;
		}
		else if (a[i] < 0.9) {
			hist[8]++;
		}
		else {
			hist[9]++;
		}
	}
	/* Print out the histgram */
	for(i=0;i<10;i++){
		printf("[%4.2f,%4.2f]:", 0.1*i, 0.1*(i+1));
		for(j=0;j<hist[i];j++) {
			// print out a symbol on every 10000th occurance
			if(j % 10000==0) {
				printf("x");
			}
		}
		printf("\n");
	}
} 
