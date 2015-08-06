ROOT = ../..

include ${ROOT}/Makefile

he:
	rm helmholtzequation; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/helmholtz; \
	make objs; \
	cd ${ROOT}/models/helmholtz-equation; \
	make helmholtzequation;

helmholtzequation:
	${CC} -o helmholtzequation HelmholtzEquation.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
