ROOT = ../..

include ${ROOT}/Makefile

pe:
	rm poissonequation; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/linsys; \
	make objs; \
	cd ${ROOT}/models/poisson-equation; \
	make poissonequation;

poissonequation:
	${CC} -o poissonequation PoissonEquation.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
