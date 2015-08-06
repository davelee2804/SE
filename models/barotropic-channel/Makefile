ROOT = ../..

include ${ROOT}/Makefile

bc:
	rm barotropicchannel; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/sw; \
	make objs; \
	cd ${ROOT}/models/barotropic-channel; \
	make barotropicchannel;

barotropicchannel:
	${CC} -o barotropicchannel BarotropicChannelSolitaryWave.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
