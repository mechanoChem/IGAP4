#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "igap4.h"
void init_linadd(double *U_, double par_dirichlet[], double par_dirichlet_temp[], void *app_)
{
    /* ================ input parameters ================ */
    App *app=(App*)app_;
    // ---- PDE
    int ndim             = app->ndim;
    int ndof             = app->ndof;
    //int *bc_type         = app->bc_type;
    //int nboun            = app->nboun;
    //int *boun            = app->boun;
    // ---- MESH
    int porder           = app->porder;
    //int *nelem           = app->nelem;
    int *nbasis          = app->nbasis;
    int *nknot           = app->nknot;
    double **knotVector  = app->knotVector;

    double dpar[ndim*ndof]; for (int i=0; i<ndim*ndof; i++) { dpar[i]=par_dirichlet[i]-par_dirichlet_temp[i]; }
    
    int ax=nbasis[1]*nbasis[2];
    int ay=nbasis[2];
    
    double *dirichlet_1d[ndim];
    double vec_p[porder];
    double mat_p[porder*porder];
    double W;
    for (int idof=0; idof<ndof; idof++)
    {
        for (int idim=0; idim<ndim; idim++)
        {
            dirichlet_1d[idim]=(double*)malloc(nbasis[idim]*sizeof(double));
            // 0th - (p-1)th deriv. at knotVector[idim][porder].
            W=knotVector[idim][porder];
            for (int i=0; i<porder; i++) { vec_p[i]=0.0; } vec_p[0]=dpar[ndim*idof+idim]*W; vec_p[1]=dpar[ndim*idof+idim];
            for (int i=0; i<porder; i++) { for (int j=0; j<porder; j++) { mat_p[porder*i+j]=evalN(knotVector[idim],nknot[idim],i,j,porder,W); }}
            matinv_dot_vec(mat_p,vec_p,porder);
            for (int i=0; i<porder; i++) { dirichlet_1d[idim][i]=vec_p[i]; }
            // 1st deriv.(=dpar[ndim*idof+idim]=const.) at knotVector[idim][i].
            for (int i=porder; i<nbasis[idim]; i++) { dirichlet_1d[idim][i]=dirichlet_1d[idim][i-1]+(knotVector[idim][i+porder]-knotVector[idim][i])/((double)porder)*dpar[ndim*idof+idim]; }
        }
        // ----
        for (int i=0; i<nbasis[0]; i++) {
            for (int j=0; j<nbasis[1]; j++) {
                for (int k=0; k<nbasis[2]; k++) {
                    U_[ndof*(ax*i+ay*j+k)+idof]+=dirichlet_1d[0][i]+dirichlet_1d[1][j]+dirichlet_1d[2][k];
        }}}
        // ----
        for (int idim=0; idim<ndim; idim++) { free(dirichlet_1d[idim]); }
    }
}

