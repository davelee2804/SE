ROOT = ../..

include ${ROOT}/Makefile

en:
	rm energetics; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/sw; \
	make objs; \
	cd ${ROOT}/tools/energetics-sw; \
	make energetics;

energetics:
	${CC} -o energetics Energetics_BTG.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
