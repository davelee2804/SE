ROOT = ../..

include ${ROOT}/Makefile

qg:
	cd ${ROOT}/src/base; \
	make objs; \
	cd ${ROOT}/src/linsys; \
	make objs; \
	cd ${ROOT}/src/qg; \
	make objs;

objs: QuasiGeostrophicEqn.o \
      TwoLayerQGEqn.o \
      libqg.a

QuasiGeostrophicEqn.o: QuasiGeostrophicEqn.cpp QuasiGeostrophicEqn.h
	${CC} -c QuasiGeostrophicEqn.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
TwoLayerQGEqn.o: TwoLayerQGEqn.cpp TwoLayerQGEqn.h
	${CC} -c TwoLayerQGEqn.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
libqg.a:
	ar -cq libqg.a *.o
	mv libqg.a ${ROOT}/lib
	cp *.h ${ROOT}/inc
	mv *.o ${ROOT}/obj
