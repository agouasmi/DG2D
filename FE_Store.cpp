#include "FE_Store.h"
#include <iostream>

FE_Store::FE_Store(Quadrature_1D *QUAD1, Quadrature_2D *QUAD2,
			FE_2D_Basis *FE_BASIS, Basis_1D* BASIS){

	/* ------ */

	FE_ORDER   = FE_BASIS->ORDER;
	FE_DIM     = FE_BASIS->DIM;

	EDGE_ORDER = BASIS->DIM-1;

	Q1_N_PTS   = QUAD1->N_PTS;
	Q2_N_PTS   = QUAD2->N_PTS;

	face_index = FE_BASIS->face_index;

	/* ------ */

	XIETA = new double[2*FE_DIM];
	
	/* ------ */

	PHI       = new double[FE_DIM*Q2_N_PTS];
	DPHI_DXI  = new double[FE_DIM*Q2_N_PTS];
	DPHI_DETA = new double[FE_DIM*Q2_N_PTS];

	/* ------ */
	
	EDGE_PHI = new double[(EDGE_ORDER+1)*Q1_N_PTS];

	/* ------ */

	Q1_WEIGHTS = QUAD1->WEIGHTS;
	Q2_WEIGHTS = QUAD2->WEIGHTS;
	Q2_POINTS  = QUAD2->POINTS;
}

void FE_Store::compute(Quadrature_1D *QUAD1, Quadrature_2D *QUAD2,
			FE_2D_Basis *FE_BASIS, Basis_1D* BASIS){

	/* ------ */
	int node;
	for (node = 0; node < FE_DIM; node++){
		FE_BASIS->xieta_at_node(&XIETA[2*node], &XIETA[2*node+1], node);	        }
	
	/* Volume terms */
	double xi, eta;
	for (int n = 0; n < Q2_N_PTS; n++){
		for (int dof = 0; dof < FE_DIM; dof++){

			xi = QUAD2->POINTS[2*n]; eta = QUAD2->POINTS[2*n+1];

			PHI[dof*Q2_N_PTS + n]       = FE_BASIS->PHI(dof,xi,eta);
			DPHI_DXI[dof*Q2_N_PTS + n]  = FE_BASIS->GRAD_PHI(dof,0,xi,eta);
			DPHI_DETA[dof*Q2_N_PTS + n] = FE_BASIS->GRAD_PHI(dof,1,xi,eta);

		}		
	}

	/* Edge terms */
	for (int n = 0; n < Q1_N_PTS; n++){
		for (int dof = 0; dof < EDGE_ORDER + 1; dof++){
	
			xi = QUAD1->POINTS[n];
			EDGE_PHI[dof*Q1_N_PTS + n] = BASIS->PHI(dof,xi);
		
		}

	}

	MASS_INV = FE_BASIS->MASS.inverse();

}

void print_vec(int* vec, int N){

	for (int i = 0; i < N; i++){
		cout << vec[i] << endl;
	}
	cout << endl;
}

void print_vec(double* vec, int N){

	for (int i = 0; i < N; i++){
		cout << vec[i] << endl;
	}
	cout << endl;
}

void FE_Store::print(){

	cout <<  "---- FE STORE ---" << endl;
	//print_vec(face_index, 4*(FE_ORDER+1));
	//print_vec(Q2_POINTS, 2*Q2_N_PTS);
	//print_vec(XIETA, 2*FE_DIM);
	cout << " PHI, DPHI_DXI, DPHI_DETA AT QUADPOINTS " << endl << endl;
	print_vec(PHI, FE_DIM*Q2_N_PTS);
	print_vec(DPHI_DXI, FE_DIM*Q2_N_PTS);
	print_vec(DPHI_DETA, FE_DIM*Q2_N_PTS);
	
	//print_vec(EDGE_PHI, (EDGE_ORDER+1)*Q1_N_PTS);

}

FE_Store::~FE_Store(){

	delete [] XIETA;
	delete [] PHI;
	delete [] DPHI_DXI;
	delete [] DPHI_DETA;
	delete [] EDGE_PHI;

	face_index = NULL;
	Q1_WEIGHTS = NULL;	
	Q2_WEIGHTS = NULL;
	Q2_POINTS  = NULL;
}
