ROOT = ../..

include ${ROOT}/Makefile

poisson:
	cd ${ROOT}/src/base; \
	make objs; \
	cd ${ROOT}/src/linsys; \
	make objs;

objs: PoissonEqn.o \
      HelmholtzEqn.o \
      StokesEqn.o \
      SchurComplement.o \
      liblinsys.a

SchurComplement.o: SchurComplement.cpp SchurComplement.h
	${CC} -c SchurComplement.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
StokesEqn.o: StokesEqn.cpp StokesEqn.h
	${CC} -c StokesEqn.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
HelmholtzEqn.o: HelmholtzEqn.cpp HelmholtzEqn.h
	${CC} -c HelmholtzEqn.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
PoissonEqn.o: PoissonEqn.cpp PoissonEqn.h
	${CC} -c PoissonEqn.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
liblinsys.a:
	ar -cq liblinsys.a *.o
	mv liblinsys.a ${ROOT}/lib
	cp *.h ${ROOT}/inc
	mv *.o ${ROOT}/obj
