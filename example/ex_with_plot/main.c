static char help[] = "gradelasttime on unitcube.\n";

#include "petscsnes.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include "igap4.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    // ---- INITIALIZE
    PetscInitialize(&argc,&argv,(char*)0,help);
    App app;




    // ---- APPLICATION
    set_ndim(&app.ndim)  ; app.nddim=1+(app.ndim)+(app.ndim)*(app.ndim+1)/2;
    set_ndof(&app.ndof)  ;
    set_torder(&app.torder);
    //                    r ,  curv,   e0 , e345 , le,   rho,  cii,        dt,    t
    double par_mat[]={ 0.25,  0.50,  500.,  250.,  0.025, 1.0,  1.0,  1.e-3*(1./2.), 0.0 }; assert_par_mat(par_mat); app.par_mat=par_mat;



    
    // ---- MESH
    mesh_uniform(&app,3,2,0,0,0);
    mesh_init(&app);




    // ---- BOUNDARY
    // F-periodic B.C.s: x=X+(F0.X+u)
    double par_periodic[] ={0.0000,0.0000,0.0000,
                            0.0000,0.0000,0.0000,
                            0.0000,0.0000,0.0000 };
    // Dirichlet/Neumann B.C.s
    //              ux,   uy,   uz, ux_n, uy_n, uz_n,       //                     ______________________
    int bc_type[]={  0,    0,    0,    1,    1,    1,       // Face 0: YZ         /                     /|
                     0,    0,    0,    1,    1,    1,       // Face 1: YZ        /          5          / |
                     0,    0,    0,    1,    1,    1,       // Face 2: ZX       /_____________________/  |
                     0,    0,    0,    1,    1,    1,       // Face 3: ZX       |                     |  |
                     0,    0,    0,    1,    1,    1,       // Face 4: XY     0 |                     | 1|    Z
                     0,    0,    0,    1,    1,    1  };    // Face 5: XY       |          2          |  /    |  Y
                                                            //                  |                     | /     | /        0: Dirichlet
                                                            //                  | ____________________|/      |/____X    1: Neumann
    double par_dirichlet[]={0.0000,0.0000,0.0000,
                            0.0000,0.0000,0.0000,
                            0.0000,0.0000,0.0000 }; // Grad u average
    double par_neumann[]  ={0.00};

    app.par_periodic = par_periodic;
    app.bc_type      = bc_type;
    app.par_dirichlet= par_dirichlet;
    app.par_neumann  = par_neumann;




    // ---- IGA
    iga_init(&app);




    // ---- ASSEMBLE/SOLVE
    int ierr;
    double *U_,*U_hist=app.U_hist;
    // I.C.s
    ierr = VecGetArray(app.U,&U_);CHKERRQ(ierr);
    appl_init(U_hist,&app);
    for (int i=0; i<app.nDof; i++) { U_hist[i]=U_hist[app.nDof+i]; }
    for (int i=0; i<app.nDof; i++) { U_[i]=U_hist[i]; }
    ierr = VecRestoreArray(app.U,&U_);CHKERRQ(ierr);
    // Time-integration
    while (par_mat[8]<5.e-4)
    {
        int converged;
        solve(&app,&converged);
        if (converged)
        {
            ierr = VecGetArray(app.U,&U_);CHKERRQ(ierr); for (int i=0; i<app.nDof; i++) { U_hist[i]=U_[i]; }
            // update solution
            par_mat[8]+=par_mat[7];
            for (int ihist=app.torder; ihist>0; ihist--) for (int i=0; i<app.nDof; i++) U_hist[app.nDof*ihist+i]=U_hist[app.nDof*(ihist-1)+i];
            ierr = VecRestoreArray(app.U,&U_);CHKERRQ(ierr);
        }
        else break;
    }




    // ---- POSTPROCESS
    compute_plot_field(&app,"plot/example.png",16,50);




    // ---- FINALIZE
    free_app(&app);
    PetscFinalize();

    return 0;
}


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          BOUNDARY CONDITION        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void eval_neumann(double *traction, double par_neumann[], int face_id, int bc_order, int idof, double xq[])
{
    int ndof=3;

    //                                          0                          |                          1                               // bc_order
    //                -------------------------------------------------------------------------------------------------------------------------
    //                        0        ,        1        ,        2        |        0        ,        1        ,        2             // idof
    //                -------------------------------------------------------------------------------------------------------------------------
    double neumann[]={        0        ,        0        ,        0        ,        0        ,        0        ,        0        ,    // face 0
                              0        ,        0        ,        0        ,        0        ,        0        ,        0        ,    // face 1
                              0        ,        0        ,        0        ,        0        ,        0        ,        0        ,    // face 2
                              0        ,        0        ,        0        ,        0        ,        0        ,        0        ,    // face 3
                              0        ,        0        ,        0        ,        0        ,        0        ,        0        ,    // face 4
                              0        ,        0        ,        0        ,        0        ,        0        ,        0         };  // face 5

    *traction=neumann[(2*ndof)*face_id+ndof*bc_order+idof];
}


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       INITIAL CONDITION/GUESS      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void appl_init(double *U_hist, void *app_)
{
    App *app=(App*)app_;
    int torder            = app->torder;
    int nDof              = app->nDof;
    for (int i=nDof; i<(torder+1)*nDof; i++) U_hist[i]=1.e-3*sin(1.0*i);
}

