#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "util.h"

int file_exist(char fname[])
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    int fexist,ierror;
    FILE *fptr;
    if (rank == 0)
    {
        fexist = ( ((fptr=fopen(fname,"rb")) != NULL) ? 1:0 );
        if (fexist == 1) { fclose(fptr); }
    }
    ierror=MPI_Bcast(&fexist,1,MPI_INT,0,MPI_COMM_WORLD);assert(ierror==MPI_SUCCESS);

    return fexist;
}



