#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "mesh.h"
#include "appctx.h"
void mesh_init(void *app_)
{
    App *app   = (App*)app_;
    int ndim   = app->ndim;
    int porder = app->porder;
    int *bc_periodic           = app->bc_periodic;
    int *nelem_global          = app->nelem_global;
    double **knotVector_global = app->knotVector_global;

    int npart[ndim],ipart[ndim];
    int *nelem=(int*) malloc(ndim*sizeof(int));
    int *nknot=(int*) malloc(ndim*sizeof(int));
    double **knotVector=(double**)malloc(ndim*sizeof(double*));
    int nsendPart,*sendPart,nrecvPart,*recvPart,ielem_displ[ndim];
    set_mpi_part(npart,ipart,nelem,ielem_displ,nknot,knotVector,&nsendPart,&sendPart,&app->sendBlock,&nrecvPart,&recvPart,&app->recvBlock,nelem_global,knotVector_global,ndim,porder,bc_periodic);

    int nboun = (int)(ipart[0]==0 && bc_periodic[0]==0)+(int)(ipart[0]==npart[0]-1 && bc_periodic[0]==0)
              + (int)(ipart[1]==0 && bc_periodic[1]==0)+(int)(ipart[1]==npart[1]-1 && bc_periodic[1]==0)
              + (int)(ipart[2]==0 && bc_periodic[2]==0)+(int)(ipart[2]==npart[2]-1 && bc_periodic[2]==0);
    int *boun=(int*)malloc(nboun*sizeof(int));
    int ia=0;
    if (ipart[0]==0 && bc_periodic[0]==0) { boun[ia++]=0; } if (ipart[0]==npart[0]-1 && bc_periodic[0]==0) { boun[ia++]=1; }
    if (ipart[1]==0 && bc_periodic[1]==0) { boun[ia++]=2; } if (ipart[1]==npart[1]-1 && bc_periodic[1]==0) { boun[ia++]=3; }
    if (ipart[2]==0 && bc_periodic[2]==0) { boun[ia++]=4; } if (ipart[2]==npart[2]-1 && bc_periodic[2]==0) { boun[ia++]=5; }

    int *nbasis=(int*) malloc(ndim*sizeof(int));
    for (int idim=0; idim<ndim; idim++) { nbasis[idim] = ( bc_periodic[idim]==0 && ipart[idim]==npart[idim]-1 ? nelem[idim]+porder : nelem[idim]); }

    set_mpi_universalIdx_temp(app,nbasis,ielem_displ,nelem_global,porder,bc_periodic);

    app->nelem        = nelem;
    app->nbasis       = nbasis;
    app->nknot        = nknot;
    app->knotVector   = knotVector;
    app->nboun        = nboun;
    app->boun         = boun;

    app->nsendPart = nsendPart;
    app->sendPart  = sendPart;
    app->nrecvPart = nrecvPart;
    app->recvPart  = recvPart;


    //free(app->bc_periodic);
    free(app->nelem_global);
    for (int idim=0; idim<ndim; idim++) { free(app->knotVector_global[idim]); }
    free(app->knotVector_global);
}
