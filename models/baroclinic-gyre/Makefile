ROOT = ../..

include ${ROOT}/Makefile

bcg:
	rm baroclinicgyre; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/sw; \
	make objs; \
	cd ${ROOT}/src/2lrl; \
	make objs; \
	cd ${ROOT}/models/baroclinic-gyre; \
	make baroclinicgyre;

baroclinicgyre:
	${CC} -o baroclinicgyre BaroclinicGyre4.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
