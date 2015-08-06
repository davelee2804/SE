ROOT = ../..

include ${ROOT}/Makefile

gr:
	rm envelope; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/tools/growth-rates; \
	make growthrates;

growthrates:
	${CC} -o growthrates GrowthRates.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
