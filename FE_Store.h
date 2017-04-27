/*   -------- Storage Class for precomputed FE values ----------- */

#ifndef FESTORE_H
#define FESTORE_H

#include <iostream>
#include <stdio.h>
#include </usr/include/eigen3/Eigen/Dense>

#include "FE_2D_basis.h"
#include "basis_1D.h"
#include "quad_1D.h"
#include "quad_2D.h"

using namespace std;
using namespace Eigen;

class FE_Store{
	
	public:
		int FE_ORDER, FE_DIM, EDGE_ORDER, Q1_N_PTS, Q2_N_PTS;

		int* face_index;

		double* XIETA;
		
		double* PHI;
		double* DPHI_DXI;
		double* DPHI_DETA;
	
		double* EDGE_PHI; 

		double* Q1_WEIGHTS;
		
		double* Q2_POINTS;
		double* Q2_WEIGHTS;

		MatrixXd MASS_INV;

		FE_Store(Quadrature_1D *QUAD1, Quadrature_2D *QUAD2, 
				FE_2D_Basis *FE_BASIS, Basis_1D *BASIS);
		void compute(Quadrature_1D *QUAD1, Quadrature_2D *QUAD2, 
				FE_2D_Basis *FE_BASIS, Basis_1D *BASIS);
		void print();
		~FE_Store();
	
};

#endif 
