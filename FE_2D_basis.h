/* ------ 2D FE Basis for DG ------ */

#ifndef FE_2D_BASIS_H
#define FE_2D_BASIS_H

#include "basis_1D.h"
#include "quad_2D.h"

#include</usr/include/eigen3/Eigen/Dense>

using namespace Eigen;

/* Base Class */
class FE_2D_Basis{
	
	public:
		int ORDER;
		int DIM;
		int NF;          

		int * face_index;                            
		MatrixXd MASS;        // Mass matrix
		
		virtual void print();
		
		FE_2D_Basis(int order, int dim, int nf);
		virtual ~FE_2D_Basis();
		
		virtual void xieta_at_node(double* xi, double* eta, unsigned int node) = 0;
    		virtual double PHI(unsigned int node, double xi, double eta) = 0;
		virtual double GRAD_PHI(unsigned int node, unsigned int dir, 
						double xi, double eta) = 0;
		
		virtual void compute_mass() = 0;
		virtual void compute_face_indices() = 0;
		virtual void setup() = 0;

};

/* Lagrange Basis on triangles */
class Lagrange_Tri: public FE_2D_Basis{

	public:
		MatrixXd C; // coefficient matrix
		double * sub_index;	
	
		void compute_C();
		void print(); 

		Lagrange_Tri(int order);
		~Lagrange_Tri(){delete sub_index;};

		int index(unsigned int r, unsigned int s);
		void xieta_at_node(double* xi, double* eta, unsigned int node);
		double PHI(unsigned int node, double xi, double eta);
		double GRAD_PHI(unsigned int node, unsigned int dir, 
						double xi, double eta);

		void compute_mass();
		void compute_face_indices();
		void setup();		

};


/* Lagrange Basis on Quad */
class Lagrange_Quad: public FE_2D_Basis{

	public:
		Basis_1D* BASIS;

		Lagrange_Quad(int order);
		~Lagrange_Quad(){delete BASIS;};

		int index(unsigned int r, unsigned int s);
		void xieta_at_node(double* xi, double* eta, unsigned int node);
		double PHI(unsigned int node, double xi, double eta);
		double GRAD_PHI(unsigned int node, unsigned int dir, 
						double xi, double eta);

		void compute_mass();
		void compute_face_indices();
		void setup();	
};

#endif
