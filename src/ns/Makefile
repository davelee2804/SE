ROOT = ../..

include ${ROOT}/Makefile

ns:
	cd ${ROOT}/src/base; \
	make objs; \
	cd ${ROOT}/src/linsys; \
	make objs;
	cd ${ROOT}/src/ns; \
	make objs;

objs: NavierStokesEqn.o \
      NeumannBC.o \
      NSEqn_KIO.o \
      libns.a

NavierStokesEqn.o: NavierStokesEqn.cpp NavierStokesEqn.h
	${CC} -c NavierStokesEqn.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
NeumannBC.o: NeumannBC.cpp NeumannBC.h
	${CC} -c NeumannBC.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
NSEqn_KIO.o: NSEqn_KIO.cpp NSEqn_KIO.h
	${CC} -c NSEqn_KIO.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
libns.a:
	ar -cq libns.a *.o
	mv libns.a ${ROOT}/lib
	cp *.h ${ROOT}/inc
	mv *.o ${ROOT}/obj
