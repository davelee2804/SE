ROOT = ../..

include ${ROOT}/Makefile

bce:
	rm baroclinicenvelope; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/linsys; \
	make objs; \
	cd ${ROOT}/src/qg; \
	make objs; \
	cd ${ROOT}/models/lee-held-91; \
	make baroclinicenvelope;

baroclinicenvelope:
	${CC} -o baroclinicenvelope BaroclinicEnvelope.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
