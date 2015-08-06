ROOT = ../..

include ${ROOT}/Makefile

st:
	rm spectralturbulence; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/ns-spectral; \
	make objs; \
	cd ${ROOT}/models/spectral-turbulence; \
	make spectralturbulence;

spectralturbulence:
	${CC} -o spectralturbulence SpectralTurbulence.cpp ${SE_OBJS} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
