#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "mesh.h"
#include "appctx.h"
void set_mpi_universalIdx_temp(void *app_, int nbasis[], int ielem_displ[], int nelem_global[], int porder, int bc_periodic[])
{
    App *app=(App*)app_;

    app->universalIdx=(int*)malloc(nbasis[0]*nbasis[1]*nbasis[2]*sizeof(int));
    int ia=0;
    for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++) {
    for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++) {
    for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++) {
        app->universalIdx[ia++]=(nelem_global[1]+((int)(bc_periodic[1]==0))*porder)*(nelem_global[2]+((int)(bc_periodic[2]==0))*porder)*(ielem_displ[0]+ibasis_x)+(nelem_global[2]+((int)(bc_periodic[2]==0))*porder)*(ielem_displ[1]+ibasis_y)+(ielem_displ[2]+ibasis_z);
    }}}
}
