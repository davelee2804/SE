ROOT = ../..

include ${ROOT}/Makefile

ns-spectral:
	cd ${ROOT}/src/base; \
	make objs; \
	cd ${ROOT}/src/qg-spectral; \
	make objs;

objs: QuasiGeostrophicSpectral.o \
      libqgs.a

QuasiGeostrophicSpectral.o: QuasiGeostrophicSpectral.cpp QuasiGeostrophicSpectral.h
	${CC} -c QuasiGeostrophicSpectral.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
libqgs.a:
	ar -cq libqgs.a *.o
	mv libqgs.a ${ROOT}/lib
	cp *.h ${ROOT}/inc
	mv *.o ${ROOT}/objs
