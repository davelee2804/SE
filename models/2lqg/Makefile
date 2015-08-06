ROOT = ../..

include ${ROOT}/Makefile

bc:
	rm baroclinicchannel; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/linsys; \
	make objs; \
	cd ${ROOT}/src/qg; \
	make objs; \
	cd ${ROOT}/models/2lqg; \
	make baroclinicchannel;

baroclinicchannel:
	${CC} -o baroclinicchannel BaroclinicChannelSolitaryWave.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
