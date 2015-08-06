ROOT = ../..

include ${ROOT}/Makefile

bg:
	rm barotropicgyre; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/sw; \
	make objs; \
	cd ${ROOT}/models/barotropic-gyre; \
	make barotropicgyre;

barotropicgyre:
	${CC} -o barotropicgyre BarotropicGyre.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
