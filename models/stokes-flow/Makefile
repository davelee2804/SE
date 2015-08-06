ROOT = ../..

include ${ROOT}/Makefile

sf:
	rm stokesflow; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/linsys; \
	make objs; \
	cd ${ROOT}/models/stokes-flow; \
	make stokesflow;

stokesflow:
	${CC} -o stokesflow StokesFlow.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
