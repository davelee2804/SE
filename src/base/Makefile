ROOT = ../..

include $(ROOT)/Makefile

objs: QuadPoint.o \
      LinAlg.o \
      Jacobi.o \
      Element.o \
      Legendre.o \
      Chebyshev.o \
      Mesh.o \
      BCs.o \
      Field.o \
      Utils.o \
      Operator.o \
      RHSOp.o \
      Advector.o \
      Vector.o \
      Matrix.o \
      SVV.o \
      libbase.a

SVV.o: SVV.cpp SVV.h Operator.o RHSOp.o Field.o BCs.o Mesh.o Element.o Jacobi.o LinAlg.o QuadPoint.o
	${CC} -c SVV.cpp ${CFLAGS}
Matrix.o: Matrix.cpp Matrix.h Vector.o Operator.o RHSOp.o Field.o BCs.o Mesh.o Element.o Jacobi.o LinAlg.o QuadPoint.o
	${CC} -c Matrix.cpp ${CFLAGS}
Vector.o: Vector.cpp Vector.h RHSOp.o Field.o BCs.o Mesh.o Element.o Jacobi.o LinAlg.o QuadPoint.o
	${CC} -c Vector.cpp ${CFLAGS}
Advector.o: Advector.cpp Advector.h Field.o BCs.o Mesh.o Element.o Jacobi.o LinAlg.o QuadPoint.o
	${CC} -c Advector.cpp ${CFLAGS}
RHSOp.o: RHSOp.cpp RHSOp.h Field.o BCs.o Mesh.o Element.o Jacobi.o LinAlg.o QuadPoint.o
	${CC} -c RHSOp.cpp ${CFLAGS}
Operator.o: Operator.cpp Operator.h Field.o BCs.o Mesh.o Element.o Jacobi.o LinAlg.o QuadPoint.o
	${CC} -c Operator.cpp ${CFLAGS}
Utils.o: Utils.cpp Utils.h Field.o BCs.o Mesh.o Legendre.o Element.o Jacobi.o LinAlg.o QuadPoint.o
	${CC} -c Utils.cpp ${CFLAGS}
Field.o: Field.cpp Field.h BCs.o Mesh.o Element.o Jacobi.o LinAlg.o QuadPoint.o
	${CC} -c Field.cpp ${CFLAGS}
BCs.o: BCs.cpp BCs.h Mesh.o Element.o Jacobi.o LinAlg.o QuadPoint.o
	${CC} -c BCs.cpp ${CFLAGS}
Mesh.o: Mesh.cpp Mesh.h Chebyshev.o Legendre.o Element.o Jacobi.o LinAlg.o QuadPoint.o
	${CC} -c Mesh.cpp ${CFLAGS}
Chebyshev.o: Chebyshev.cpp Chebyshev.h Element.o Jacobi.o LinAlg.o QuadPoint.o
	${CC} -c Chebyshev.cpp ${CFLAGS}
Legendre.o: Legendre.cpp Legendre.h Element.o Jacobi.o LinAlg.o QuadPoint.o
	${CC} -c Legendre.cpp ${CFLAGS}
Element.o: Element.cpp Element.h Jacobi.o LinAlg.o QuadPoint.o
	${CC} -c Element.cpp ${CFLAGS}
Jacobi.o: Jacobi.cpp Jacobi.h
	${CC} -c Jacobi.cpp ${CFLAGS}
LinAlg.o: LinAlg.cpp LinAlg.h
	${CC} -c LinAlg.cpp ${CFLAGS}
QuadPoint.o: QuadPoint.cpp QuadPoint.h
	${CC} -c QuadPoint.cpp ${CFLAGS}

libbase.a:
	ar -cq libbase.a *.o
	mv libbase.a ${ROOT}/lib
	cp *.h ${ROOT}/inc
	mv *.o ${ROOT}/obj
