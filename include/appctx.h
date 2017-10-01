#if !defined (APPCTX)
#define APPCTX

#include "petscsnes.h"
/* ---------------- declaration (structure) ---------------- */
typedef struct{
    // ---- BVP
    int ndim;
    int nddim;
    int ndof;
    int torder;
    double *par_mat;
    int *bc_type;
    double *par_periodic;
    double *par_dirichlet;
    double *par_neumann;
    int nboun;
    int *boun;
    // ---- MESH
    int porder;
    int *nelem;
    int *nbasis;
    int *nknot;
    double **knotVector;
    // ---- IGA
    Vec U;
    Vec Residual;
    Mat Tangent;
    SNES snes;
    int nquad;
    double *cquad;
    double *wquad;
    int nselfIdx;
    int *selfIdx;
    int nsendPart;
    int *sendPart;
    int nsendIdx;
    int *sendIdx;
    int *sendPtr;
    int nrecvPart;
    int *recvPart;
    int nrecvIdx;
    int *recvIdx;
    int *recvPtr;
    int *globalIdx;
    double *U_hist;

    // ---- data storage
    int nDof;
    int *nDof_array;
    int *iDof_displ_array;
    int nDof_global;
    int *universalIdx;

    // MESH -> IGA (transient storage)
    int *sendBlock;
    int *recvBlock;

    // Global MESH Info (transient storage)
    int *bc_periodic;
    int *nelem_global;
    double **knotVector_global;

} App;

void free_app(void *app_);

#endif
