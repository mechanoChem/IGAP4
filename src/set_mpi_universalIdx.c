#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "igap4.h"
void set_mpi_universalIdx(void *app_)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);
    
    App *app=(App*)app_;
    int ndof=app->ndof;

    int *universalIdx_self=(int*)malloc(app->nselfIdx*sizeof(int));
    int ia=0, ja=0;
    for (int ibasis_x=0; ibasis_x<app->nbasis[0]; ibasis_x++) {
    for (int ibasis_y=0; ibasis_y<app->nbasis[1]; ibasis_y++) {
    for (int ibasis_z=0; ibasis_z<app->nbasis[2]; ibasis_z++) {
        for (int idof=0; idof<ndof; idof++) { universalIdx_self[ia++]=ndof*(app->universalIdx[ja])+idof; } ja++;
    }}}
    free(app->universalIdx);
    if (rank==0) { app->universalIdx=(int*)malloc(app->nDof_global*sizeof(int)); } else { app->universalIdx=(int*)malloc(0); }
    
    assert(MPI_Gatherv(universalIdx_self,app->nDof,MPI_INT,app->universalIdx,app->nDof_array,app->iDof_displ_array,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    free(universalIdx_self);
}
