ROOT = ../..

include ${ROOT}/Makefile

qg:
	rm testqg; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/linsys; \
	make objs; \
	cd ${ROOT}/src/qg; \
	make objs; \
	cd ${ROOT}/models/test-qg; \
	make testqg;

testqg:
	${CC} -o testqg TestQG.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
