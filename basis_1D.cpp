/* ---- Function definitions for class 1D_basis ---- */

#include "basis_1D.h"
#include <cstddef>
#include <iostream>

using namespace std;

/* Lagrange Basis methods */
Lagrange::Lagrange(int order): Basis_1D(order){

	POINTS = new double[DIM];
	if (DIM > 1){
		double dx = 2./ (double) (DIM-1);
		for (int i = 0; i < DIM; i++){
			POINTS[i] = -1. + (double) i * dx;
		}
	}else
	{
		POINTS[0] = 0.;
	}

}

void Lagrange::print(){

	cout << endl << "---- 1D BASIS ----" << endl << endl;
	cout << "Ref Points for the 1D Basis:" << endl;
	for (int i = 0; i < DIM; i++){
		cout << POINTS[i] << ", ";
	}
	cout << endl << endl;

}

double Lagrange::PHI(unsigned int node, double xi){
	
	int n; double out = 1.;
	double xi_ref = POINTS[node];
	double xi_temp;

	for (n = 0; n < DIM; n++){
		if (n != node){
			xi_temp = POINTS[n];
			out *= (xi - xi_temp) / (xi_ref - xi_temp);
		}
	}
		
	return out;
}

double Lagrange::GRAD_PHI(unsigned int node, double xi){
        
        int n,p;
        double denom = 1.; double temp = 1.; double out = 0.;

        double xi_ref = POINTS[node];
        double xi_temp;

        // Compute the common denominator first
        for (n = 0; n < DIM; n++){
                if (n != node){
                        xi_temp = POINTS[n];
                        denom *= 1 / (xi_ref - xi_temp);
                }
        }

        // Do the rest  
        for (p = 0; p < DIM; p++){
                if (p != node){
                        for (n = 0; n < DIM; n++){
                                if ((n != node)&&(n != p)){
                                        xi_temp = POINTS[n];
                                        temp *= (xi - xi_temp);
                                }
                        }
                out += temp;
                temp = 1.;
                }
        }

        return out*denom;
}



/* Legendre Basis methods */
Legendre::Legendre(int order): Basis_1D(order){

	C = MatrixXd::Zero(DIM,DIM);

	if (order < 1){
		C(0,0) = 1.;
	}
	
	if (order < 2){
		C(0,1) = 0.; 
		C(1,1) = 1.;
	}

	if (order < 3){
		C(0,2) = -1./2.; 
		C(1,2) = 0.;
		C(2,2) = 3./2.;
	}

	if (order < 4){
		C(0,3) = 0.;
		C(1,3) = -3./2.;
		C(2,3) = 0.;
		C(3,3) = 5./2.;
	}

	if (order < 5){
		C(0,4) = 3./8.;
		C(1,4) = 0.;
		C(2,4) = -30./8.;
		C(3,4) = 0.;
		C(4,4) = 35./8.;
	}

	if (order < 6){
		C(0,5) = 0.;
		C(1,5) = 15./8.;
		C(2,5) = 0.;
		C(3,5) = -70./8.;
		C(4,5) = 0.;
		C(5,5) = 63./8.;
	}

	if (order < 7){
		C(0,6) = -5./16.;
		C(1,6) = 0.;
		C(2,6) = 105./16.;
		C(3,6) = 0.;
		C(4,6) = -315./16.;
		C(5,6) = 0.;
		C(6,6) = 231./16.;
	}

	if (order < 8){
		C(0,7) = 0.;
		C(1,7) = -35./16.;
		C(2,7) = 0.;
		C(3,7) = -315./16.;
		C(4,7) = 0.;
		C(5,7) = -693./16.;
		C(6,7) = 0.;
		C(7,7) = 429./16.;
	}
	
}

double Legendre::PHI(unsigned int mode, double xi){

	double out = 0.;
	for (int k = 0; k < DIM; k++){
		out += C(k,mode)*pow(xi,k);
	}
	return out;
}

double Legendre::GRAD_PHI(unsigned int mode, double xi){

	double out = 0.;
	for (int k = 1; k < DIM; k++){
		out += double(k) * C(k,mode)*pow(xi,k-1);
	} 

}

void Legendre::print(){

	cout << C << endl;

}
