#if !defined (BOUNDARY)
#define BOUNDARY

void eval_neumann(double *traction, double par_neumann[], int face_id, int bc_order, int idof, double xq[]);
void appl_dirichlet(double *dirichlet, double par[], int face_id, int bc_order, int idof, double *knotVector_[], int nknot_[], int nbasis_[], int ndim, int porder);
void appl_init(double *U_hist, void *app_);

#endif
