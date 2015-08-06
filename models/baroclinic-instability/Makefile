ROOT = ../..

include ${ROOT}/Makefile

bci:
	rm baroclinicenvelope; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/linsys; \
	make objs; \
	cd ${ROOT}/src/qg; \
	make objs; \
	cd ${ROOT}/models/baroclinic-instability; \
	make baroclinicinstability;

baroclinicinstability:
	${CC} -o barocliniceinstability BaroclinicInstability.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
