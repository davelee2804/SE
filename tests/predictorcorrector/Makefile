ROOT = ../..

include ${ROOT}/Makefile

pc:
	rm predictorcorrector; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/sw; \
        make objs; \
	cd ${ROOT}/src/2lrl; \
        make objs; \
	cd ${ROOT}/tests/predictorcorrector; \
	make predictorcorrector;

predictorcorrector:
	${CC} -o predictorcorrector TestPredictorCorrector.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
