ROOT = ../..

include ${ROOT}/Makefile

pv:
	rm potentialvorticitykelvinhelmholtz; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/poisson; \
        make objs; \
	cd ${ROOT}/src/pv; \
        make objs; \
	cd ${ROOT}/models/potential-vorticity; \
	make potentialvorticitykelvinhelmholtz;

potentialvorticitykelvinhelmholtz:
	${CC} -o potentialvorticitykelvinhelmholtz PotentialVorticityKelvinHelmholtz.cpp ${SE_OBJS} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
