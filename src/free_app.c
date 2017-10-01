#include "petscsnes.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "appctx.h"
void free_app(void *app_)
{
    App *app = (App*)app_;
    
    free(app->nelem);
    free(app->nbasis);

    for (int idim=0; idim<app->ndim; idim++) { free(app->knotVector[idim]); } // allocated in "set_mpi_part.c"
    free(app->sendPart);                                                 //
    free(app->recvPart);                                                 //
    free(app->selfIdx);                                                  // allocated in "set_mpi_comm.c"
    free(app->sendIdx);                                                  //
    free(app->sendPtr);                                                  //
    free(app->recvIdx);                                                  //
    free(app->recvPtr);                                                  //
    free(app->cquad);
    free(app->wquad);
    free(app->globalIdx);
    free(app->nDof_array);
    free(app->iDof_displ_array);
    free(app->universalIdx);

    free(app->U_hist);
    
    VecDestroy(&app->U);
    VecDestroy(&app->Residual);
    MatDestroy(&app->Tangent);
    SNESDestroy(&app->snes);
}
