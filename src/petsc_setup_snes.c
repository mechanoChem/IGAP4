#include "petscsnes.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "igap4.h"

void petsc_setup_snes(void *app_)
{
    App            *app = (App*)app_;
    KSP            ksp;
    PC             pc;
    SNESLineSearch ls;
    // ---- SNES
    SNESCreate(MPI_COMM_WORLD,&app->snes);
    SNESSetTolerances(app->snes,1.e-10,1.e-20,1.e-20,1,1e9);
    SNESSetFunction(app->snes,app->Residual,compute_Residual,app);
    SNESSetJacobian(app->snes,app->Tangent,app->Tangent,compute_Tangent,app);
    SNESSetType(app->snes,SNESNEWTONLS);
    SNESGetLineSearch(app->snes,&ls);
    SNESLineSearchSetType(ls,SNESLINESEARCHBT); //SNESLINESEARCHBT, SNESLINESEARCHL2, SNESLINESEARCHCP, SNESLINESEARCHBASIC, SNESLINESEARCHSHELL
    SNESLineSearchSetOrder(ls,SNES_LINESEARCH_ORDER_CUBIC); //LINEAR, QUADRATIC, CUBIC
    // ---- KSP/PC
    SNESGetKSP(app->snes,&ksp); 
    KSPGetPC(ksp,&pc);
    KSPSetInitialGuessNonzero(ksp,1);
    KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,1e8,1e8);
//if DYNAM == 0
//    KSPSetType(ksp,KSPMINRES); // The operator and the preconditioner must be symmetric and the preconditioner must be positive definite for this method.
//    KSPSetPCSide(ksp,PC_LEFT); //
//    PCSetType(pc,PCJACOBI);
//    PCJacobiSetUseAbs(pc,1);
//elif DYNAM == 1
    KSPSetType(ksp,KSPGMRES);
    KSPGMRESSetRestart(ksp,10000); // eats memory. 2e4 doesn't work on edison; srun error.
    KSPGMRESSetOrthogonalization(ksp,KSPGMRESModifiedGramSchmidtOrthogonalization);//KSPGMRESClassicalGramSchmidtOrthogonalization
    KSPSetPCSide(ksp,PC_RIGHT);//PC_RIGHT, PC_LEFT, PC_SYMMETRIC
    PCSetType(pc,PCASM);
    PCASMSetType(pc,PC_ASM_RESTRICT);//PC_ASM_BASIC, PC_ASM_INTERPOLATE, PC_ASM_RESTRICT, PC_ASM_NONE
/*    KSPSetType(ksp,KSPCG);
    KSPCGSetType(ksp,KSP_CG_SYMMETRIC);
    PCSetType(pc,PCJACOBI);
    PCJacobiSetUseAbs(pc,1);*/
//endif
    // ---- overwrite with runtime option
    SNESSetFromOptions(app->snes);
}
