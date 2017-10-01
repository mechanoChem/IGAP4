#include "petscsnes.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "igap4.h"

void petsc_setup_arrays(void *app_)
{
    App *app=(App *)app_;
    int nDof        = app->nDof;
    int nDof_global = app->nDof_global;
    
    // ---- Vec
    VecCreateMPI(MPI_COMM_WORLD,nDof,nDof_global,&app->Residual);
    VecDuplicate(app->Residual,&app->U);
    // ---- Mat
    int *d_nz=(int*) malloc(nDof*sizeof(int)),*o_nz=(int*) malloc(nDof*sizeof(int));
    set_nz(d_nz,o_nz,app->nbasis,app->ndof,app->porder,app->nsendPart,app->sendBlock,app->nrecvPart,app->recvBlock);
    MatCreateAIJ(MPI_COMM_WORLD,nDof,nDof,nDof_global,nDof_global,0,d_nz,0,o_nz,&app->Tangent);
    free(d_nz);free(o_nz);
    MatSetOption(app->Tangent,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);
    MatSetOption(app->Tangent,MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE);
    MatSetOption(app->Tangent,MAT_SUBSET_OFF_PROC_ENTRIES,PETSC_TRUE);
//if DYNAM == 0
//    MatSetOption(Tangent,MAT_SYMMETRIC,PETSC_TRUE);
//    MatSetOption(Tangent,MAT_SYMMETRY_ETERNAL,PETSC_TRUE);
//endif
}
