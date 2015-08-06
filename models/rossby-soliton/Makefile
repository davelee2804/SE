ROOT = ../..

include ${ROOT}/Makefile

rs:
	rm rossbysoliton; \
	cd ${ROOT}/src/base; \
        make objs; \
	cd ${ROOT}/src/sw; \
	make objs; \
	cd ${ROOT}/models/rossby-soliton; \
	make rossbysoliton;

rossbysoliton:
	${CC} -o rossbysoliton RossbySoliton.cpp ${SE_OBJ} \
        ${PETSC_LIB} ${HDF5_LIB} ${CFLAGS} ${SE_LIB} ${SE_INC}
