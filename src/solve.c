#include "petscsnes.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include "igap4.h"
void solve(void *app_, int *converged)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    App *app=(App*)app_;
    clock_t time_start, time_end;

    KSP ksp; SNESGetKSP(app->snes,&ksp);
    int                 ksp_niter;
    SNESConvergedReason snes_reason;

    // ---- subiterate
    double Fnorm_new,Fnorm_old=999.;
    for (int isubiter=0; isubiter<10000; isubiter++)
    {
        time_start=clock();
        SNESSolve(app->snes,NULL,app->U);
        KSPGetIterationNumber(ksp,&ksp_niter);
        SNESGetConvergedReason(app->snes,&snes_reason);
        SNESGetFunctionNorm(app->snes,&Fnorm_new); if ( double_abs(Fnorm_old-Fnorm_new)/Fnorm_old < 1.e-6 ) { snes_reason=SNES_DIVERGED_LINE_SEARCH; } Fnorm_old=Fnorm_new;
        time_end=clock();
        if (rank==0) { printf("%6d: Fnorm = %e, #lin.ite.= %6d, (%08.2f[min])\n",isubiter+1,Fnorm_new,ksp_niter,((float)(time_end-time_start))/CLOCKS_PER_SEC/60.0); }
        if (snes_reason != SNES_DIVERGED_MAX_IT) { break; }
    }
    if (snes_reason<=0 && rank==0) { printf("%s.\n",SNESConvergedReasons[snes_reason]); }
    if (snes_reason>0) { *converged=1; } else { *converged=0; }
}
