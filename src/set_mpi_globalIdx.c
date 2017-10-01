#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "igap4.h"
void set_mpi_globalIdx(int *globalIdx, int *iDof_displ_array, int nselfIdx, int *selfIdx, int nsendPart, int *sendPart, int nsendIdx, int *sendIdx, int *sendPtr, int nrecvPart, int *recvPart, int nrecvIdx, int *recvIdx, int *recvPtr, int nbasis[], int ndof)
{
    MPI_Request request_send[nsendPart],request_recv[nrecvPart];
    MPI_Status status;
    int ierror;
    int ia=0,ibasis; 
    int iDof_displ;
    ierror=MPI_Scatter(iDof_displ_array,1,MPI_INT,&iDof_displ,1,MPI_INT,0,MPI_COMM_WORLD);assert(ierror==MPI_SUCCESS);
    int *globalIdx_self = (int*) malloc(nselfIdx*sizeof(int));
    for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
    for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
    for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
        ibasis=(nbasis[1]*nbasis[2])*ibasis_x+nbasis[2]*ibasis_y+ibasis_z;
        for (int idof=0; idof<ndof; idof++) { globalIdx_self[ia++]=iDof_displ+ndof*ibasis+idof; }
    }}}
    for (int i=0; i<nselfIdx; i++) { globalIdx[selfIdx[i]]=globalIdx_self[i]; }
    free(globalIdx_self);
    int *globalIdx_send=(int*) malloc(nsendIdx*sizeof(int));
    int *globalIdx_recv=(int*) malloc(nrecvIdx*sizeof(int));
    for (int i=0; i<nsendIdx; i++) { globalIdx_send[i]=globalIdx[sendIdx[i]]; }
    for (int isendPart=0; isendPart<nsendPart; isendPart++) { ierror=MPI_Isend(globalIdx_send+sendPtr[isendPart],sendPtr[isendPart+1]-sendPtr[isendPart],MPI_INT,sendPart[isendPart],1,MPI_COMM_WORLD,request_send+isendPart);assert(ierror==MPI_SUCCESS); }
    for (int irecvPart=0; irecvPart<nrecvPart; irecvPart++) { ierror=MPI_Irecv(globalIdx_recv+recvPtr[irecvPart],recvPtr[irecvPart+1]-recvPtr[irecvPart],MPI_INT,recvPart[irecvPart],1,MPI_COMM_WORLD,request_recv+irecvPart);assert(ierror==MPI_SUCCESS); }
    for (int isendPart=0; isendPart<nsendPart; isendPart++) assert(MPI_Wait(request_send+isendPart,&status)==MPI_SUCCESS);
    for (int irecvPart=0; irecvPart<nrecvPart; irecvPart++) assert(MPI_Wait(request_recv+irecvPart,&status)==MPI_SUCCESS);
    for (int i=0; i<nrecvIdx; i++) { globalIdx[recvIdx[i]]=globalIdx_recv[i]; }
    free(globalIdx_send);
    free(globalIdx_recv);
}
