#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "igap4.h"
void set_nz(int *d_nz, int *o_nz, int nbasis[], int ndof, int porder, int nsendPart, int *sendBlock, int nrecvPart, int *recvBlock)
{
    int sd[]={0,0,0}; for (int i=0; i<nsendPart; i++) { if (sendBlock[i]<3) sd[sendBlock[i]]=1; }
    int rv[]={0,0,0}; for (int i=0; i<nrecvPart; i++) { if (recvBlock[i]<3) rv[recvBlock[i]]=1; }
    int ia=0;
    for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++) {
    for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++) {
    for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++) {
    for (int idof=0; idof<ndof; idof++) {
        d_nz[ia]=ndof*( ( ibasis_x < porder ? ibasis_x : porder )+1+( ibasis_x > (nbasis[0]-1)-porder ? (nbasis[0]-1)-ibasis_x:porder) )
                     *( ( ibasis_y < porder ? ibasis_y : porder )+1+( ibasis_y > (nbasis[1]-1)-porder ? (nbasis[1]-1)-ibasis_y:porder) )
                     *( ( ibasis_z < porder ? ibasis_z : porder )+1+( ibasis_z > (nbasis[2]-1)-porder ? (nbasis[2]-1)-ibasis_z:porder) );
        o_nz[ia]=ndof*( ( (sd[2]==0 && ibasis_x < porder) ? ibasis_x : porder )+1+( (rv[2]==0 && ibasis_x > (nbasis[0]-1)-porder ) ? (nbasis[0]-1)-ibasis_x : porder ) )
                     *( ( (sd[1]==0 && ibasis_y < porder) ? ibasis_y : porder )+1+( (rv[1]==0 && ibasis_y > (nbasis[1]-1)-porder ) ? (nbasis[1]-1)-ibasis_y : porder ) )
                     *( ( (sd[0]==0 && ibasis_z < porder) ? ibasis_z : porder )+1+( (rv[0]==0 && ibasis_z > (nbasis[2]-1)-porder ) ? (nbasis[2]-1)-ibasis_z : porder ) )
                 -d_nz[ia];
        ia++;
    }}}}
}

// if send to/recv from self: d_nz=d_nz+o_nz, o_nz=0.



