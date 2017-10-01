#if !defined (DOMAIN)
#define DOMAIN

void set_ndim(int *ndim);
void set_ndof(int *ndof);
void set_torder(int *torder);

void eval_residual(double *residual, double *u, double *par);
void eval_tangent(double *tangent, double *u, double *par);


void set_nfield(int *nfield);
void eval_field(double *field, double *u, double *par);


void assert_par_mat(double *par);
#endif
