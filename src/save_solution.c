#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "igap4.h"
void save_solution(double *U_, char fname[], void *app_)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    App *app=(App*)app_;
    int nDof              = app->nDof;
    int *nDof_array       = app->nDof_array;
    int *iDof_displ_array = app->iDof_displ_array;
    int nDof_global       = app->nDof_global;
    int *universalIdx     = app->universalIdx;

    double *U_universal     ; if (rank==0) { U_universal     =(double*) malloc(nDof_global*sizeof(double)); } else { U_universal     =(double*) malloc(0); }
    double *U_universal_temp; if (rank==0) { U_universal_temp=(double*) malloc(nDof_global*sizeof(double)); } else { U_universal_temp=(double*) malloc(0); }
    assert(MPI_Gatherv(U_,nDof,MPI_DOUBLE,U_universal_temp,nDof_array,iDof_displ_array,MPI_DOUBLE,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    if (rank==0)
    {
        for (int i=0; i<nDof_global; i++) { U_universal[universalIdx[i]]=U_universal_temp[i]; }
        FILE *fptr=fopen(fname,"wb");fwrite(U_universal,sizeof(double),nDof_global,fptr);fclose(fptr);
    }

    free(U_universal);
    free(U_universal_temp);
}
