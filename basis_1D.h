/* ------- 1D Polynomial basis [-1 1] ------- */

#ifndef BASIS_1D_H
#define BASIS_1D_H

#include </usr/include/eigen3/Eigen/Dense>

using namespace Eigen;

class Basis_1D{
	
	public:
		int DIM;                
		
		virtual double PHI(unsigned int node, double xi) = 0;
		virtual double GRAD_PHI(unsigned int node, double xi) = 0;
		virtual void print() = 0;
		
		Basis_1D(int order){DIM = order+1;};
		virtual ~Basis_1D(){};

};

class Lagrange: public Basis_1D{

	public:
		double * POINTS;

		double PHI(unsigned int node, double xi);
		double GRAD_PHI(unsigned int node, double xi);
		void print();

		Lagrange(int order);
		~Lagrange(){delete POINTS;};
};

class Legendre: public Basis_1D{

	public: 
		MatrixXd C;

		double PHI(unsigned int mode, double xi);
		double GRAD_PHI(unsigned int mode, double xi);
		void print();		

		Legendre(int order);
		
};

#endif
