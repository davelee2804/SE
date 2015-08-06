ROOT = ../..

include ${ROOT}/Makefile

2lrl:
	cd ${ROOT}/src/base; \
	make objs; \
	cd ${ROOT}/src/sw; \
	make objs; 
	cd ${ROOT}/src/2lrl; \
	make objs;

objs: TwoLayerSW_Energetics.o \
      TwoLayerSW.o \
      BarotropicEqn.o \
      SW2L.o \
      2LRLOps.o \
      lib2lrl.a

TwoLayerSW_Energetics.o: TwoLayerSW_Energetics.cpp TwoLayerSW_Energetics.h TwoLayerSW.o BarotropicEqn.o 2LRLOps.o
	${CC} -c TwoLayerSW_Energetics.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
TwoLayerSW.o: TwoLayerSW.cpp TwoLayerSW.h 2LRLOps.o
	${CC} -c TwoLayerSW.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
BarotropicEqn.o: BarotropicEqn.cpp BarotropicEqn.h 2LRLOps.o
	${CC} -c BarotropicEqn.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
SW2L.o: SW2L.cpp SW2L.h 2LRLOps.o
	${CC} -c SW2L.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
2LRLOps.o: 2LRLOps.cpp 2LRLOps.h
	${CC} -c 2LRLOps.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
lib2lrl.a:
	ar -cq lib2lrl.a *.o
	mv lib2lrl.a ${ROOT}/lib
	cp *.h ${ROOT}/inc
	mv *.o ${ROOT}/obj
