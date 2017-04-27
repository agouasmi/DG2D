/* ---- Function definitions for class 2D_FE_basis ---- */

#include "FE_2D_basis.h"
#include "math.h"

#include <iostream>

using namespace std;
using namespace Eigen;

/* Base Class */
FE_2D_Basis::FE_2D_Basis(int order, int dim, int nf){

	ORDER = order;
	DIM = dim;
	NF = nf;

	face_index = new int[NF*(ORDER + 1)];
	MASS = MatrixXd::Zero(DIM,DIM);
	
}

FE_2D_Basis::~FE_2D_Basis(){
	
	delete [] face_index;
		
}

void FE_2D_Basis::print(){

	cout << endl << "---- FE BASIS ----" << endl << endl;
	cout << "FE Basis has been loaded, the indices of the face nodes are: " << endl;

	for (int face = 0; face < NF; face++){
		cout << "[ " ;
		for (int i = 0; i < (ORDER + 1); i++){
			cout << face_index[face*(ORDER + 1) + i] << " " ;
		}
		cout << "] ";
	}
	cout << endl << endl;
	cout << "The Mass matrix is given by: " << endl << MASS << endl;
	cout << "There are " << DIM << " DOFS " << endl;
}

/* Lagrange basis on triangles */

int get_dim(int order){

	return (int) (order+1)*(order+2)/2;
}

Lagrange_Tri::Lagrange_Tri(int order)
                     :FE_2D_Basis(order, get_dim(order), 3){

	C = MatrixXd::Zero(DIM,DIM);
	
}

void Lagrange_Tri::compute_C(){

	int r, s, node, R, S, k;
	double xi, eta;

	MatrixXd A(DIM,DIM);

	for (s = 0; s < ORDER + 1; s++){
		for (r = 0; r < ORDER + 1 - s; r++){
			
			node = index(r,s);
			xi   = (double) r/ (double) ORDER;
			eta  = (double) s/ (double) ORDER;

			for (S = 0; S < ORDER + 1; S++){
				for (R = 0; R < ORDER + 1 - S; R++){
				
					k         = index(R,S);		
					A(node,k) = pow(xi,R)*pow(eta,S);

				}
			}

		}
	}
	C = A.inverse();

}


void Lagrange_Tri::print(){

	FE_2D_Basis::print();
	cout << "the coefficient matrix is given by:" << endl << C << endl;

}

int Lagrange_Tri::index(unsigned int r, unsigned int s){

	int k = r;
	for (int s_ = 0; s_ < s; s_++){
		k += ORDER + 1 - s_;
	}

	return k;
}

void Lagrange_Tri::xieta_at_node(double *xi, double *eta, unsigned int node){

	int r,s;
	
	r = sub_index[2*node];
	s = sub_index[2*node+1];

	if (ORDER != 0){
		*xi  = (double) r/ORDER;
		*eta = (double) s/ORDER;
	}
	else {
		*xi  = 0.33;
		*eta = 0.33;
	}
}

double Lagrange_Tri::PHI(unsigned int node, double xi, double eta){

	double out = 0.;
	int r, s, k;
	
	for (s = 0; s < ORDER + 1; s++){
		for (r = 0; r < ORDER + 1 - s; r++){
			
			k = index(r,s);
			out += C(k,node) * pow(xi,r) * pow(eta,s);
		}
	}
		
	return out;

}

double Lagrange_Tri::GRAD_PHI(unsigned int node, unsigned int dir, double xi, double eta){
	
	double out = 0.;
	int r, s, k;
	
	if (dir == 0){ // Diff wrt xi
		for (s = 0; s < ORDER + 1; s++){
			for (r = 1; r < ORDER + 1 - s; r++){

				k = index(r,s);
				out += r * C(k,node) * pow(xi,r-1) * pow(eta,s);

			}
		}
	}	 
	else // Diff wrt eta
	{  
		for (s = 1; s < ORDER + 1; s++){
			for (r = 0; r < ORDER + 1 - s; r++){
			
				k = index(r,s);
				out += s * C(k,node) * pow(xi,r) * pow(eta,s-1);

			}
		}

	}
	return  out;

}



void Lagrange_Tri::compute_mass(){
	
	int i, j, n;
	double xi, eta;

	Dunavant Quad(2*ORDER);
	
	for (i = 0; i < DIM; i++){
		for (j = 0; j < DIM; j++){
			
			MASS(i,j) = 0.;
			for (n = 0; n < Quad.N_PTS; n++){

				xi  = Quad.POINTS[2*n];
				eta = Quad.POINTS[2*n+1];
				MASS(i,j) += Quad.WEIGHTS[n]*PHI(i,xi,eta)*PHI(j,xi,eta);

			}			

		}	
	}

}

void Lagrange_Tri::compute_face_indices(){

	int r, s, i;

	// First face
	r = 0; s = 0;
	for (i = 0; i < ORDER + 1; i++){
		face_index[0*(ORDER + 1) + i] = index(r,s);
		r++;
	}

	// Second face
	r = ORDER; s = 0;
	for (i = 0; i < ORDER + 1; i++){
		face_index[1*(ORDER + 1) + i] = index(r,s);
		r--;
		s++;
	}

	// Third face	
	r = 0; s = ORDER;
	for (i = 0; i < ORDER + 1; i++){
		face_index[2*(ORDER + 1) + i] = index(r,s);
		s--;
	}

}

void Lagrange_Tri::setup(){

	sub_index = new double[2*DIM];
	int r,s;

	for (s = 0; s < ORDER + 1; s++){
		for (r = 0; r < ORDER + 1 - s; r++){

			sub_index[2*index(r,s)]   = r;
			sub_index[2*index(r,s)+1] = s;
		}
	}	

	compute_C();
	compute_mass();
	compute_face_indices();

}


/* Lagrange basis on Quads (Tensor) */

Lagrange_Quad::Lagrange_Quad(int order): FE_2D_Basis(order, (order+1)*(order+1), 4){
	BASIS = new Lagrange(order);
}

int Lagrange_Quad::index(unsigned int r, unsigned int s){

	return s*(ORDER+1) + r;
}

void Lagrange_Quad::xieta_at_node(double *xi, double *eta, unsigned int node){

	int r, s;

	r = node % (ORDER+1);
	s = (int) (node - r)/(ORDER + 1);
	
	if (ORDER != 0){
		*xi  = -1 + 2*(double) r/ORDER;
		*eta = -1 + 2*(double) s/ORDER;
	}
	else {
		*xi  = 0;
		*eta = 0;
	}
}

double Lagrange_Quad::PHI(unsigned int node, double xi, double eta){

	double out = 0.;
	int r, s;

	r = node % (ORDER+1);
	s = (int) (node - r)/(ORDER + 1);	

	out = BASIS->PHI(r,xi) * BASIS->PHI(s,eta);
			
	return out;

}

double Lagrange_Quad::GRAD_PHI(unsigned int node, unsigned int dir, double xi, double eta){
	
	double out = 0.;
	int r, s;
	
	r = node % (ORDER+1);
	s = (int) (node - r)/(ORDER + 1);
	
	if (dir == 0){ // Diff wrt xi
		out += BASIS->GRAD_PHI(r,xi) * BASIS->PHI(s,eta);
	}	 
	else // Diff wrt eta
	{  
		out += BASIS->PHI(r,xi) * BASIS->GRAD_PHI(s,eta);
	}

	return  out;
}

void Lagrange_Quad::compute_mass(){
	
	int i, j, n;
	double xi, eta;

	Gauss quad1(ORDER+1);
	Tensor Quad(&quad1);
	
	for (i = 0; i < DIM; i++){
		for (j = 0; j < DIM; j++){
			
			MASS(i,j) = 0.;
			for (n = 0; n < Quad.N_PTS; n++){

				xi  = Quad.POINTS[2*n];
				eta = Quad.POINTS[2*n+1];
				MASS(i,j) += Quad.WEIGHTS[n]*PHI(i,xi,eta)*PHI(j,xi,eta);

			}			

		}	
	}

}

void Lagrange_Quad::compute_face_indices(){

	int r, s, i;

	// First face
	r = 0; s = 0;
	for (i = 0; i < ORDER + 1; i++){
		face_index[0*(ORDER + 1) + i] = index(r,s);
		r++;
	}

	// Second face
	r = ORDER; s = 0;
	for (i = 0; i < ORDER + 1; i++){
		face_index[1*(ORDER + 1) + i] = index(r,s);
		s++;
	}

	// Third face
	r = ORDER; s = ORDER;
	for (i = 0; i < ORDER + 1; i++){
		face_index[2*(ORDER + 1) + i] = index(r,s);
		r--;
	}

	// Fourth face	
	r = 0; s = ORDER;
	for (i = 0; i < ORDER + 1; i++){
		face_index[3*(ORDER + 1) + i] = index(r,s);
		s--;
	}

}

void Lagrange_Quad::setup(){

	compute_mass();
	compute_face_indices();

}

