#include "mpi.h"
#include <stdio.h>

int main (int argc, char **argv)
{
    int myid, numprocs;
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &myid);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    printf ("Hello World! Process %d of %d\n", myid, numprocs);
    MPI_Finalize ();
    return 0;
}

