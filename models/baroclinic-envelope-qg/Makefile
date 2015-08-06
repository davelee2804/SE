ROOT = ../..

include ${ROOT}/Makefile

bce:
	rm baroclinicenvelope; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/sw; \
	make objs; \
	cd ${ROOT}/src/qg; \
	make objs; \
	cd ${ROOT}/models/baroclinic-envelope-qg; \
	make baroclinicenvelope;

baroclinicenvelope:
	${CC} -o baroclinicenvelope BaroclinicEnvelope.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
