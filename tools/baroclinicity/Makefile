ROOT = ../..

include ${ROOT}/Makefile

bc:
	rm baroclinicity; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/tools/baroclinicity; \
	make baroclinicity;

baroclinicity:
	${CC} -o baroclinicity Baroclinicity.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
