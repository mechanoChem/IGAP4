# IGAP4
Isogeometric Analysis Program for (initial) boundary value problems that involve higher-order spatial derivatives. 

## Overview
IGAP4 is written in C. It requires PETSc [https://www.mcs.anl.gov/petsc/] for linear/nonlinear solvers. It can optionally use MathGL [http://mathgl.sourceforge.net/doc_en/Main.html] for quick plottings. 
IGAP4 can be used to numerically solve (initial) boundary value problems that involve second-order spatial derivatives in the weak forms; examples include (transient) gradient elasticity and the Cahn-Hilliard equation. 

## Installation on Linux

1) Install PETSc 3.7.x [https://www.mcs.anl.gov/petsc/]. Set environmental variables PETSC_DIR and PETSC_ARCH as required by PETSc.

2) Set environmental variable, IGAP4_DIR, on the command line:
> export IGAP4_DIR=/path/to/your/igap4/dir

3) Compile the source and make shared libraries:
> cd ${IGAP4_DIR}; make

4) If MathGL is installed and one desires to use the quick plotting functionality of IGAP4, do:
> cd ${IGAP4_DIR}/src/mgl; make

## Tests

Compile the example application, "gradelasttimets", for testing (this will take a few minutes):
> cd ${IGAP4_DIR}/application/gradelasttimets; make

Example without MathGL plottings (with 8 processes):
> cd ${IGAP4_DIR}/example/ex_wo_plot; make; mpiexec -n 8 ./main

Example with MathGL plottings (with 8 processes):
> cd ${IGAP4_DIR}/example/ex_with_plot; make; mpiexec -n 8 ./main

## Contributors

This code has been developed at the Computational Physics Group at the University of Michigan [http://www.umich.edu/~compphys/index.html].

- Koki Sagiyama (Lead Developer)
- Krishna Garikipati

## Acknowledgements

The development of this code has been supported by the following:

- Dept of Energy (DoE and labs) : Software Center for Predictive Theory and Modeling (DE-SC0008637) 
- National Science Foundation : Integrated Computational Framework for Designing Dynamically Controlled Alloy-Oxide Heterostructures (DMR1436154)

## License

GNU LESSER GENERAL PUBLIC LICENSE v3.0
