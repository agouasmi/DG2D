EXEC = DG2D

OBJECTS = main.o DG_2D.o Element.o Geometry.o Equation.o FE_2D_basis.o basis_1D.o quad_1D.o quad_2D.o FE_Store.o

SOURCES = main.cpp DG_2D.cpp Element.cpp Gemoetry.cpp Equation.cpp FE_2D_basis.cpp basis_1D.cpp quad_1D.cpp quad_2D.cpp FE_Store.cpp

HEADERS = DG_2D.h Element.h Geometry.h Equation.h FE_2D_basis.h basis_1D.h quad_1D.h quad_2D.h utils.h anim.h FE_Store.h

CFLAGS = -c -std=c++11 -O3 -I/usr/include #-I/usr/dislin #-g

LFLAGS = #-L/usr/local/dislin -ldiscpp -larmadillo -ldislin -lXt -lm 

CC = g++


${EXEC}: ${OBJECTS}
	${CC} ${OBJECTS} ${LFLAGS} -o ${EXEC}

main.o: main.cpp DG_2D.h 
	${CC} ${CFLAGS} main.cpp DG_2D.h  

DG_2D.o: DG_2D.cpp DG_2D.h 
	${CC} ${CFLAGS} DG_2D.cpp DG_2D.h  

Element.o: Element.cpp Element.h
	${CC} ${CFLAGS} Element.cpp Element.h

Equation.o: Equation.cpp Equation.h
	${CC} ${CFLAGS} Equation.cpp Equation.h

quad_1D.o: quad_1D.cpp quad_1D.h
	${CC} ${CFLAGS} quad_1D.cpp quad_1D.h

quad_2D.o: quad_2D.cpp quad_2D.h
	${CC} ${CFLAGS} quad_2D.cpp quad_2D.h

FE_2D_basis.o: FE_2D_basis.cpp FE_2D_basis.h
	${CC} ${CFLAGS} FE_2D_basis.cpp FE_2D_basis.h

basis_1D.o: basis_1D.cpp basis_1D.h
	${CC} ${CFLAGS} basis_1D.cpp basis_1D.h 

FE_Store.o: FE_Store.cpp FE_Store.h
	${CC} ${CFLAGS} FE_Store.cpp FE_Store.h

Geometry.o: Geometry.cpp Geometry.h
	${CC} ${CFLAGS} Geometry.cpp Geometry.h


clean:
	rm -r ${OBJECTS} ${EXEC}
	rm *gch*
	rm Output/solution*

run:
	./${EXEC}

