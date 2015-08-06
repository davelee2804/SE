ROOT = ../..

include ${ROOT}/Makefile

ns-spectral:
	cd ${ROOT}/src/base; \
	make objs; \
	cd ${ROOT}/src/ns-spectral; \
	make objs;

objs: NavierStokesSpectral.o \
      libnss.a

NavierStokesSpectral.o: NavierStokesSpectral.cpp NavierStokesSpectral.h
	${CC} -c NavierStokesSpectral.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
libnss.a:
	ar -cq libnss.a *.o
	mv libnss.a ${ROOT}/lib
	cp *.h ${ROOT}/inc
	mv *.o ${ROOT}/objs
