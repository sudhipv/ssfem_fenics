LLAPACK = -llapack 
LBLAS = -lblas
MAIN = ./a.out
COMP = mpif90
COMP_OPS = -O2

PETSC_DIR=/Users/ajit/Downloads/SOFTWAREs/petsc-3.7.6
PETSC_ARCH=arch-darwin-c-opt
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

LIST = common.o myCommon.o variationalform.o PETScAssembly.o assembly.o PETScommon.o PETScSolvers.o main.o

all : $(LIST) chkopts
	-${FLINKER} -o $(MAIN) $(LIST) -lpetsc 
	${RM} *.o *.mod
