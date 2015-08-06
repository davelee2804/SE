ROOT = ../..

include ${ROOT}/Makefile

rk:
	rm rossbykelvin irrotational; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/sw; \
	make objs; \
	cd ${ROOT}/models/rossby-kelvin; \
	make rossbykelvin;

ir:
	rm rossbykelvin irrotational; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/sw; \
	make objs; \
	cd ${ROOT}/models/rossby-kelvin; \
	make irrotational;

rossbykelvin:
	${CC} -o rossbykelvin RossbyKelvin.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}

irrotational:
	${CC} -o irrotational Irrotational.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
