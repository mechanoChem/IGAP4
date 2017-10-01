#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "mesh.h"
void set_mpi_part(
int npart[], 
int ipart[], 
int nelem[], 
int ielem_displ[], 
int nknot[], 
double *knotVector[], 
int *nsendPart_,
int **sendPart_,
int **sendBlock_,
int *nrecvPart_,
int **recvPart_,
int **recvBlock_,
int nelem_global[], 
double *knotVector_global[], 
int ndim, 
int porder, 
int bc_periodic[])
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    int nepp_approx;
    nepp_approx=round(pow((double)(nelem_global[0]*nelem_global[1]*nelem_global[2])/(double)nproc,1./3.));
    for (int i=nelem_global[0]/nepp_approx+1; i>0; i--)  { if (nproc%i ==0) { npart[0]=i; break;} } 
    nepp_approx=round(pow((double)(nelem_global[1]*nelem_global[2])/(double)(nproc/npart[0]),1./2.));
    for (int i=nelem_global[1]/nepp_approx+1; i>0; i--)  { if (nproc/npart[0]%i ==0) { npart[1]=i; break;} }
    npart[2]=nproc/npart[0]/npart[1];
    //for (int idim=0; idim<ndim; idim++) { assert(nelem_global[idim]/npart[idim]>=2*porder); }
    if (rank==0) { printf("npart[]={%d,%d,%d}\n",npart[0],npart[1],npart[2]); }

    ipart[0]=rank/(npart[1]*npart[2]);
    ipart[1]=(rank%(npart[1]*npart[2]))/npart[2];
    ipart[2]=rank%npart[2];

    int *nelem_w_array;       if (rank==0) { nelem_w_array      =(int*) malloc(nproc*sizeof(int)); } else { nelem_w_array      =(int*) malloc(0); }
    int *nknot_w_array;       if (rank==0) { nknot_w_array      =(int*) malloc(nproc*sizeof(int)); } else { nknot_w_array      =(int*) malloc(0); }
    int *ielem_w_displ_array; if (rank==0) { ielem_w_displ_array=(int*) malloc(nproc*sizeof(int)); } else { ielem_w_displ_array=(int*) malloc(0); }
    for (int idim=0; idim<ndim; idim++)
    {
        if (rank==0)
        {
            int *nelem_w_1d;
            int *ielem_w_displ_1d;
            int ipart_w;
            nelem_w_1d      =(int*) malloc(npart[idim]*sizeof(int));
            ielem_w_displ_1d=(int*) malloc(npart[idim]*sizeof(int));
            for (int i=0; i<npart[idim]; i++) { nelem_w_1d[i]=nelem_global[idim]/npart[idim]; }
            for (int i=0; i<nelem_global[idim]%npart[idim]; i++) { nelem_w_1d[i]++; }
            ielem_w_displ_1d[0]=0;
            for (int i=1; i<npart[idim]; i++) { ielem_w_displ_1d[i]=ielem_w_displ_1d[i-1]+nelem_w_1d[i-1]; }
            for (int iproc=0; iproc<nproc; iproc++)
            {
                switch (idim) {
                case 0: ipart_w=iproc/(npart[1]*npart[2]); break;
                case 1: ipart_w=(iproc%(npart[1]*npart[2]))/npart[2]; break;
                case 2: ipart_w=iproc%npart[2]; break;
                default: exit(EXIT_FAILURE);
                }
                nelem_w_array[iproc]=nelem_w_1d[ipart_w];
                nknot_w_array[iproc]=nelem_w_1d[ipart_w]+1+2*porder;
                ielem_w_displ_array[iproc]=ielem_w_displ_1d[ipart_w];
            }
            free(nelem_w_1d);
            free(ielem_w_displ_1d);
        }
        assert(MPI_Scatter(nelem_w_array,1,MPI_INT,nelem+idim,1,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
        assert(MPI_Scatter(nknot_w_array,1,MPI_INT,nknot+idim,1,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
        assert(MPI_Scatter(ielem_w_displ_array,1,MPI_INT,ielem_displ+idim,1,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
        knotVector[idim]=(double*) malloc(nknot[idim]*sizeof(double));
        assert(MPI_Scatterv(knotVector_global[idim],nknot_w_array,ielem_w_displ_array,MPI_DOUBLE,knotVector[idim],nknot[idim],MPI_DOUBLE,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    }
    free(nelem_w_array);
    free(nknot_w_array);
    free(ielem_w_displ_array);

    int ia,cond0,cond1,cond2;
    // ---- send
    cond0=( bc_periodic[0]==1 || ipart[0] != 0 );
    cond1=( bc_periodic[1]==1 || ipart[1] != 0 );
    cond2=( bc_periodic[2]==1 || ipart[2] != 0 );
    int nsendPart=(int)(                   cond2 )
                 +(int)(          cond1          )
                 +(int)( cond0                   )
                 +(int)(          cond1 && cond2 )
                 +(int)( cond0          && cond2 )
                 +(int)( cond0 && cond1          )
                 +(int)( cond0 && cond1 && cond2 );
    int *sendPart =(int*)malloc(nsendPart*sizeof(int));
    int *sendBlock=(int*)malloc(nsendPart*sizeof(int));
    ia=0;
    if (                   cond2 ) { sendPart[ia]=npart[1]*npart[2]*((ipart[0]  +npart[0])%npart[0])+npart[2]*((ipart[1]  +npart[1])%npart[1])+((ipart[2]-1+npart[2])%npart[2]); sendBlock[ia]=0; ia++; }
    if (          cond1          ) { sendPart[ia]=npart[1]*npart[2]*((ipart[0]  +npart[0])%npart[0])+npart[2]*((ipart[1]-1+npart[1])%npart[1])+((ipart[2]  +npart[2])%npart[2]); sendBlock[ia]=1; ia++; }
    if ( cond0                   ) { sendPart[ia]=npart[1]*npart[2]*((ipart[0]-1+npart[0])%npart[0])+npart[2]*((ipart[1]  +npart[1])%npart[1])+((ipart[2]  +npart[2])%npart[2]); sendBlock[ia]=2; ia++; }
    if (          cond1 && cond2 ) { sendPart[ia]=npart[1]*npart[2]*((ipart[0]  +npart[0])%npart[0])+npart[2]*((ipart[1]-1+npart[1])%npart[1])+((ipart[2]-1+npart[2])%npart[2]); sendBlock[ia]=3; ia++; }
    if ( cond0          && cond2 ) { sendPart[ia]=npart[1]*npart[2]*((ipart[0]-1+npart[0])%npart[0])+npart[2]*((ipart[1]  +npart[1])%npart[1])+((ipart[2]-1+npart[2])%npart[2]); sendBlock[ia]=4; ia++; }
    if ( cond0 && cond1          ) { sendPart[ia]=npart[1]*npart[2]*((ipart[0]-1+npart[0])%npart[0])+npart[2]*((ipart[1]-1+npart[1])%npart[1])+((ipart[2]  +npart[2])%npart[2]); sendBlock[ia]=5; ia++; }
    if ( cond0 && cond1 && cond2 ) { sendPart[ia]=npart[1]*npart[2]*((ipart[0]-1+npart[0])%npart[0])+npart[2]*((ipart[1]-1+npart[1])%npart[1])+((ipart[2]-1+npart[2])%npart[2]); sendBlock[ia]=6; ia++; }
    // ---- recv
    cond0=( bc_periodic[0]==1 || ipart[0] != npart[0]-1 );
    cond1=( bc_periodic[1]==1 || ipart[1] != npart[1]-1 );
    cond2=( bc_periodic[2]==1 || ipart[2] != npart[2]-1 );
    int nrecvPart=(int)(                   cond2 )
                 +(int)(          cond1          )
                 +(int)( cond0                   )
                 +(int)(          cond1 && cond2 )
                 +(int)( cond0          && cond2 )
                 +(int)( cond0 && cond1          )
                 +(int)( cond0 && cond1 && cond2 );
    int *recvPart =(int*)malloc(nrecvPart*sizeof(int));
    int *recvBlock=(int*)malloc(nrecvPart*sizeof(int));
    ia=0;
    if (                   cond2 ) { recvPart[ia]=npart[1]*npart[2]*((ipart[0]           )%npart[0])+npart[2]*((ipart[1]           )%npart[1])+((ipart[2]+1         )%npart[2]); recvBlock[ia]=0; ia++; }
    if (          cond1          ) { recvPart[ia]=npart[1]*npart[2]*((ipart[0]           )%npart[0])+npart[2]*((ipart[1]+1         )%npart[1])+((ipart[2]           )%npart[2]); recvBlock[ia]=1; ia++; }
    if ( cond0                   ) { recvPart[ia]=npart[1]*npart[2]*((ipart[0]+1         )%npart[0])+npart[2]*((ipart[1]           )%npart[1])+((ipart[2]           )%npart[2]); recvBlock[ia]=2; ia++; }
    if (          cond1 && cond2 ) { recvPart[ia]=npart[1]*npart[2]*((ipart[0]           )%npart[0])+npart[2]*((ipart[1]+1         )%npart[1])+((ipart[2]+1         )%npart[2]); recvBlock[ia]=3; ia++; }
    if ( cond0          && cond2 ) { recvPart[ia]=npart[1]*npart[2]*((ipart[0]+1         )%npart[0])+npart[2]*((ipart[1]           )%npart[1])+((ipart[2]+1         )%npart[2]); recvBlock[ia]=4; ia++; }
    if ( cond0 && cond1          ) { recvPart[ia]=npart[1]*npart[2]*((ipart[0]+1         )%npart[0])+npart[2]*((ipart[1]+1         )%npart[1])+((ipart[2]           )%npart[2]); recvBlock[ia]=5; ia++; }
    if ( cond0 && cond1 && cond2 ) { recvPart[ia]=npart[1]*npart[2]*((ipart[0]+1         )%npart[0])+npart[2]*((ipart[1]+1         )%npart[1])+((ipart[2]+1         )%npart[2]); recvBlock[ia]=6; ia++; }
    
    // ----
    *nsendPart_ = nsendPart;
    *sendPart_  = sendPart;
    *sendBlock_ = sendBlock;
    *nrecvPart_ = nrecvPart;
    *recvPart_  = recvPart;
    *recvBlock_ = recvBlock;
}
