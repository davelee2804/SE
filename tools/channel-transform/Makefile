ROOT = ../..

include ${ROOT}/Makefile

ct:
	rm channeltransform; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/transform; \
	make objs; \
	cd ${ROOT}/tools/channel-transform; \
	make channeltransform;

channeltransform:
	${CC} -o channeltransform CT.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
