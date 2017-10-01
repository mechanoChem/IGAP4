#if !defined (MESH)
#define MESH

void mesh_init(void *app_);
void set_mpi_part(int npart[], int ipart[], int nelem[], int ielem_displ[], int nknot[], double *knotVector[], int *nsendPart_, int **sendPart_, int **sendBlock_, int *nrecvPart_, int **recvPart_, int **recvBlock_, int nelem_global[], double *knotVector_global[], int ndim, int porder, int bc_periodic[]);
void mesh_uniform(void *app_, int mref, int porder, int bc_periodic_x, int bc_periodic_y, int bc_periodic_z);
void eval_Bspline(double N[], double *kV[], int p, double xi[]);
double evalN(double *kV, int nk, int k, int i, int p, double xi);
void set_mpi_universalIdx_temp(void *app_, int nbasis[], int ielem_displ[], int nelem_global[], int porder, int bc_periodic[]);

#endif
