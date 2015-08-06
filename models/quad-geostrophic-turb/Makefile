ROOT = ../..

include ${ROOT}/Makefile

qgt:
	rm quasigeostrophicturb; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/qg-spectral; \
	make objs; \
	cd ${ROOT}/models/quasi-geostrophic-turb; \
	make quasigeostrophicturb;

quasigeostrophicturb:
	${CC} -o quasigeostrophicturb QuasiGeostrophicTurb.cpp ${SE_OBJS} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
