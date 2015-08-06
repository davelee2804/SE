ROOT = ../..

include ${ROOT}/Makefile

sw:
	cd ${ROOT}/src/base; \
	make objs; \
	cd ${ROOT}/src/sw; \
	make objs;

objs: Energetics_SW.o \
      ShallowWaterEqn.o \
      LinearShallowWaterEqn.o \
      SWOps.o \
      libsw.a

Energetics_SW.o: Energetics_SW.cpp Energetics_SW.h ShallowWaterEqn.o
	${CC} -c Energetics_SW.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
ShallowWaterEqn.o: ShallowWaterEqn.cpp ShallowWaterEqn.h
	${CC} -c ShallowWaterEqn.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
LinearShallowWaterEqn.o: LinearShallowWaterEqn.cpp LinearShallowWaterEqn.h
	${CC} -c LinearShallowWaterEqn.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
SWOps.o: SWOps.cpp SWOps.h
	${CC} -c SWOps.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
libsw.a:
	ar -cq libsw.a *.o
	mv libsw.a ${ROOT}/lib
	cp *.h ${ROOT}/inc
	mv *.o ${ROOT}/obj
