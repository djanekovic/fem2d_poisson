PETSC_DIR = /home/darko/FEM/petsc

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

.DEFAULT_GOAL := poisson

poisson: poisson.o chkopts
	-${CLINKER} poisson.o -g ${PETSC_LIB}
	${RM} poisson.o

run: a.out square.msh
	${PETSC_DIR}/${PETSC_ARCH}/bin/mpiexec -n 2 ./a.out -dim 2

mesh: square.geo
	gmsh -2 square.geo
