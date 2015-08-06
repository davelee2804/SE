ROOT = ../..

include ${ROOT}/Makefile

tv:
	rm taylorvortex; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/linsys; \
	make objs; \
	cd ${ROOT}/src/ns; \
	make objs; \
	cd ${ROOT}/models/taylor-vortex; \
	make taylorvortex;

taylorvortex:
	${CC} -o taylorvortex TaylorVortex.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
