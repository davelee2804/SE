ROOT = ../..

include ${ROOT}/Makefile

kh:
	rm kelvinhelmholtz; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/linsys; \
	make objs; \
	cd ${ROOT}/src/ns; \
	make objs; \
	cd ${ROOT}/models/kelvin-helmholtz; \
	make kelvinhelmholtz;

kelvinhelmholtz:
	${CC} -o kelvinhelmholtz KelvinHelmholtz.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
