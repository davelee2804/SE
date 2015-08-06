ROOT = ../..

include ${ROOT}/Makefile

svv:
	rm testSVV; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/linsys; \
	make objs;
	cd ${ROOT}/tests/testSVV; \
	make testSVV;

testSVV:
	${CC} -o testSVV TestSVV.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
