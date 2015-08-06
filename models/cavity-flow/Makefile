ROOT = ../..

include ${ROOT}/Makefile

cf:
	rm cavityflow; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/linsys; \
	make objs; \
	cd ${ROOT}/src/ns; \
	make objs; \
	cd ${ROOT}/models/cavity-flow; \
	make cavityflow;

cavityflow:
	${CC} -o cavityflow CavityFlow.cpp ${SE_OBJS} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
