#include "petscsnes.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "igap4.h"
void iga_init(void* app_)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    App *app=(App*)app_;
    int ndof      = app->ndof;
    int torder    = app->torder;
    int porder    = app->porder;
    int *nelem    = app->nelem;
    int *nbasis   = app->nbasis;
    int nsendPart = app->nsendPart;
    int *sendPart = app->sendPart;
    int nrecvPart = app->nrecvPart;
    int *recvPart = app->recvPart;

    // ----
    int nselfIdx,*selfIdx,nsendIdx,*sendIdx,*sendPtr,nrecvIdx,*recvIdx,*recvPtr;
    set_mpi_comm(&nselfIdx,&selfIdx,&nsendIdx,&sendIdx,&sendPtr,&nrecvIdx,&recvIdx,&recvPtr,nsendPart,app->sendBlock,nrecvPart,app->recvBlock,ndof,porder,nelem,nbasis);
    app->nselfIdx = nselfIdx;
    app->selfIdx  = selfIdx;
    app->nsendIdx = nsendIdx;
    app->sendIdx  = sendIdx;
    app->sendPtr  = sendPtr;
    app->nrecvIdx = nrecvIdx;
    app->recvIdx  = recvIdx;
    app->recvPtr  = recvPtr;

    // ----
    int nDof=ndof*(nbasis[0]*nbasis[1]*nbasis[2]);
    int *nDof_array;          if (rank==0) { nDof_array      =(int*) malloc(nproc*sizeof(int)); }     else { nDof_array      =(int*) malloc(0); }
    int *iDof_displ_array;    if (rank==0) { iDof_displ_array=(int*) malloc((nproc+1)*sizeof(int)); } else { iDof_displ_array=(int*) malloc(0); }
    assert(MPI_Gather(&nDof,1,MPI_INT,nDof_array,1,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    if (rank==0) { iDof_displ_array[0]=0; for (int iproc=1; iproc<nproc+1; iproc++) { iDof_displ_array[iproc]=iDof_displ_array[iproc-1]+nDof_array[iproc-1]; } }
    int nDof_global;
    if (rank==0) { nDof_global=iDof_displ_array[nproc]; } assert(MPI_Bcast(&nDof_global,1,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    app->nDof             = nDof;
    app->nDof_array       = nDof_array;
    app->iDof_displ_array = iDof_displ_array;
    app->nDof_global      = nDof_global;

    // ----
    petsc_setup_arrays(app); free(app->sendBlock); free(app->recvBlock);
    petsc_setup_snes(app);
    app->U_hist=(double*)malloc((torder+1)*nDof*sizeof(double));

    // ---- globalIdx (maps local idx to global idx) required by Petsc
    int *globalIdx      = (int*) malloc((nselfIdx+nrecvIdx)*sizeof(int));
    set_mpi_globalIdx(globalIdx,iDof_displ_array,nselfIdx,selfIdx,nsendPart,sendPart,nsendIdx,sendIdx,sendPtr,nrecvPart,recvPart,nrecvIdx,recvIdx,recvPtr,nbasis,ndof);
    app->globalIdx        = globalIdx;

    // ---- universalIdx (for partition independant data storage)
    set_mpi_universalIdx(app);

    // ----
    int     nquad=4;
    double *cquad=(double*)malloc(nquad*sizeof(double));
    double *wquad=(double*)malloc(nquad*sizeof(double));
    legendre_handle(cquad,wquad,nquad,0.,1.);
    app->nquad     = nquad;
    app->cquad     = cquad;
    app->wquad     = wquad;	    
}
