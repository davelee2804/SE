ROOT = ../..

include ${ROOT}/Makefile

pv:
	cd ${ROOT}/src/base; \
	make objs; \
	cd ${ROOT}/src/pv; \
	make objs;

objs: PotentialVorticityEqn.o \
      libpv.a

PotentialVorticityEqn.o: PotentialVorticityEqn.cpp PotentialVorticityEqn.h
	${CC} -c PotentialVorticityEqn.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
libpv.a:
	ar -cq libpv.a *.o
	mv libpv.a ${ROOT}/lib
	cp *.h ${ROOT}/inc
	mv *.o ${ROOT}/objs
