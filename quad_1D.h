#ifndef QUAD_1D_H
#define QUAD_1D_H

#include "math.h"

/* Class for 1D quadrature rules on [-1, 1] */
class Quadrature_1D{

	public:

		int ORDER;
		int N_PTS;
	
		double * WEIGHTS;
		double * POINTS; 
		
		Quadrature_1D(int npts);
		void print(); 
		
		virtual ~Quadrature_1D();
		virtual void set_quad() = 0;		

};


class Gauss: public Quadrature_1D{
	
	public:
		void set_quad();		
		Gauss(int npts);
		
};

#endif
