ROOT = ../..

include ${ROOT}/Makefile

transform:
	cd ${ROOT}/src/base; \
	make objs; \
	cd ${ROOT}/src/transform; \
	make objs;

objs: FourierChebyshev.o \
      Transform.o \
      ChebTrans.o \
      libtransform.a

FourierChebyshev.o: FourierChebyshev.cpp FourierChebyshev.h ChebTrans.o
	${CC} -c FourierChebyshev.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
Transform.o: Transform.cpp Transform.h
	${CC} -c Transform.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
ChebTrans.o: ChebTrans.cpp ChebTrans.h
	${CC} -c ChebTrans.cpp ${CFLAGS} ${SE_LIB} ${SE_INC}
libtransform.a:
	ar -cq libtransform.a *.o
	mv libtransform.a ${ROOT}/lib
	cp *.h ${ROOT}/inc
	mv *.o ${ROOT}/obj
