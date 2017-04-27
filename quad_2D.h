#ifndef QUAD_2D_H
#define QUAD_2D_H

#include "quad_1D.h"
#include "math.h"

/* Class for 2D quadrature rule */
class Quadrature_2D{

	public:
		int ORDER;
		int N_PTS;
	
		double * WEIGHTS;
		double * POINTS;
 
		Quadrature_2D(int order, int npts); 
		virtual ~Quadrature_2D();

		void print();
};

/* Triangles */
class Dunavant: public Quadrature_2D{

	public:
		Dunavant(int order);
		int get_npts(int order);	
};

/* Quadrilaterals */
class Tensor: public Quadrature_2D{

	public:
		Tensor(Quadrature_1D *q);

}; 

#endif
