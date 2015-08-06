ROOT = ../..

include ${ROOT}/Makefile

bte:
	rm barotropicenvelope; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/sw; \
	make objs; \
	cd ${ROOT}/models/barotropic-envelope-sw; \
	make barotropicenvelope;

barotropicenvelope:
	${CC} -o barotropicenvelope BarotropicEnvelope.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
