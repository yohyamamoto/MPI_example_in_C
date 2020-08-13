## MPI example in C

This program uses a random number generator and the Box-Muller (polar form) transformation to generate normally distributed random numbers.
I used MPI to demonstrate some for-loops can get a benefit of parallel computings (although this is probably not the best place to implement it).

Compile this with
  mpicc sample_code_mpi.c
and run it with 
  mpirun ./a.out -np 4 
or 
  mpirun -hostfile machinefile -np 4 ./a.out

© March 19th, 2014	Yoh Yamamoto
