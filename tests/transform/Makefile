ROOT = ../..

include ${ROOT}/Makefile

tr:
	rm transform; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/transform; \
        make objs; \
	cd ${ROOT}/tests/transform; \
	make transform;

transform:
	${CC} -o transform TestTransform.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
