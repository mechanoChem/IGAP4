#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "igap4.h"
void set_mpi_comm(
int *nselfIdx_,
int **selfIdx_,
int *nsendIdx_,
int **sendIdx_,
int **sendPtr_,
int *nrecvIdx_,
int **recvIdx_,
int **recvPtr_,
int nsendPart,
int *sendBlock,
int nrecvPart,
int *recvBlock,
int ndof,
int porder,
int nelem[],
int nbasis[]
)
{
    int nbasis_y_active=nelem[1]+porder;
    int nbasis_z_active=nelem[2]+porder;
    
    // ---- 

    int ibasis,ia;
    // ---- self
    int nselfIdx=ndof*(nbasis[0]*nbasis[1]*nbasis[2]);
    int *selfIdx=(int*)malloc(nselfIdx*sizeof(int));
    ia=0;
    for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
    for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
    for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
        ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
        for (int idof=0; idof<ndof; idof++) { selfIdx[ia++]=ndof*ibasis+idof; }
    }}}

    // ---- send
    int nsendIdx=0;
    int *sendPtr=(int*)malloc((nsendPart+1)*sizeof(int)); sendPtr[0]=nsendIdx;
    for (int i=0; i<nsendPart; i++)
    {
        switch (sendBlock[i]) {
        case 0: nsendIdx+=ndof*(nbasis[0]*nbasis[1]*porder); break;
        case 1: nsendIdx+=ndof*(nbasis[0]*nbasis[2]*porder); break;
        case 2: nsendIdx+=ndof*(nbasis[1]*nbasis[2]*porder); break;
        case 3: nsendIdx+=ndof*(nbasis[0]*porder*porder); break;
        case 4: nsendIdx+=ndof*(nbasis[1]*porder*porder); break;
        case 5: nsendIdx+=ndof*(nbasis[2]*porder*porder); break;
        case 6: nsendIdx+=ndof*(porder*porder*porder); break;
        }
        sendPtr[i+1]=nsendIdx;
    }
    int *sendIdx=(int*)malloc(nsendIdx*sizeof(int));
    ia=0;
    for (int i=0; i<nsendPart; i++)
    {
        switch (sendBlock[i]) {
        case 0: // (face0)
            for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 1: // (face1)
            for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<porder  ; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 2: // (face2)
            for (int ibasis_x=0; ibasis_x<porder  ; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 3: // (edge0)
            for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<porder  ; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 4: // (edge1)
            for (int ibasis_x=0; ibasis_x<porder  ; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 5: // (edge2)
            for (int ibasis_x=0; ibasis_x<porder  ; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<porder  ; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 6: // (corner)
            for (int ibasis_x=0; ibasis_x<porder  ; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<porder  ; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        }
    }

    // ---- recv
    int nrecvIdx=0;
    int *recvPtr=(int*)malloc((nrecvPart+1)*sizeof(int)); recvPtr[0]=nrecvIdx;
    for (int i=0; i<nrecvPart; i++)
    {
        switch (recvBlock[i]) {
        case 0: nrecvIdx+=ndof*(nbasis[0]*nbasis[1]*porder); break;
        case 1: nrecvIdx+=ndof*(nbasis[0]*nbasis[2]*porder); break;
        case 2: nrecvIdx+=ndof*(nbasis[1]*nbasis[2]*porder); break;
        case 3: nrecvIdx+=ndof*(nbasis[0]*porder*porder); break;
        case 4: nrecvIdx+=ndof*(nbasis[1]*porder*porder); break;
        case 5: nrecvIdx+=ndof*(nbasis[2]*porder*porder); break;
        case 6: nrecvIdx+=ndof*(porder*porder*porder); break;
        }
        recvPtr[i+1]=nrecvIdx;
    }
    int *recvIdx=(int*)malloc(nrecvIdx*sizeof(int));
    ia=0;
    for (int i=0; i<nrecvPart; i++)
    {
        switch (recvBlock[i]) {
        case 0: // (face0)
            for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
            for (int ibasis_z=nbasis[2]; ibasis_z<nbasis[2]+porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 1: // (face1)
            for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
            for (int ibasis_y=nbasis[1]; ibasis_y<nbasis[1]+porder  ; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 2: // (face2)
            for (int ibasis_x=nbasis[0]; ibasis_x<nbasis[0]+porder  ; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 3: // (edge0)
            for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
            for (int ibasis_y=nbasis[1]; ibasis_y<nbasis[1]+porder  ; ibasis_y++){
            for (int ibasis_z=nbasis[2]; ibasis_z<nbasis[2]+porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 4: // (edge1)
            for (int ibasis_x=nbasis[0]; ibasis_x<nbasis[0]+porder  ; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
            for (int ibasis_z=nbasis[2]; ibasis_z<nbasis[2]+porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 5: // (edge2)
            for (int ibasis_x=nbasis[0]; ibasis_x<nbasis[0]+porder  ; ibasis_x++){
            for (int ibasis_y=nbasis[1]; ibasis_y<nbasis[1]+porder  ; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 6: // (corner)
            for (int ibasis_x=nbasis[0]; ibasis_x<nbasis[0]+porder  ; ibasis_x++){
            for (int ibasis_y=nbasis[1]; ibasis_y<nbasis[1]+porder  ; ibasis_y++){
            for (int ibasis_z=nbasis[2]; ibasis_z<nbasis[2]+porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        }
    }
    // ----
    *nselfIdx_      = nselfIdx;
    *nsendIdx_      = nsendIdx;
    *nrecvIdx_      = nrecvIdx;

    *selfIdx_       = selfIdx;
    *sendIdx_       = sendIdx;
    *sendPtr_       = sendPtr;
    *recvIdx_       = recvIdx;
    *recvPtr_       = recvPtr;
}
