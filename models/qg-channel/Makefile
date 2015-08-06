ROOT = ../..

include ${ROOT}/Makefile

bc:
	rm barotropicchannel; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/linsys; \
	make objs; \
	cd ${ROOT}/src/qg; \
	make objs; \
	cd ${ROOT}/models/qg-channel; \
	make barotropicchannel;

barotropicchannel:
	${CC} -o barotropicchannel BarotropicChannel.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
