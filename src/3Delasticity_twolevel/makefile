MAIN = ./a.out

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

LIST = common.o myCommon.o variationalform.o PETScAssembly.o assembly.o PETScommon.o PETScSolvers.o main.o

all : $(LIST) chkopts
	-${FLINKER} -o $(MAIN) $(LIST) ${PETSC_KSP_LIB}

# 	${RM} -f *.o *.mod
#
