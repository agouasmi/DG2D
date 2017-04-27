/*   -------- 2D Finite Element Class ----------- */

#ifndef ELEMENT_H
#define ELEMENT_H

#include <iostream>
#include <stdio.h>
#include </usr/include/eigen3/Eigen/Dense>

#include "FE_Store.h"
#include "Equation.h"

using namespace std;
using namespace Eigen;

class Element{
	
	public:
		int M, ORDER, DOFS_NUMBER, NF; 
        	
		// Geometry info
		double             size;
		MatrixXd              N;
		VectorXd        X, Y, L;
	
		// Inverse of Mass matrix in global space
		MatrixXd        IMASS;

		// Faces average max wave speeds + CFL number
		//double *smax;
		//double CFL;

		// Data structures
		double *dofs_backup;           // Temporary var for low storage RK4
		double *dofs;                  // DOF vector 
		double *dofs_rhs;
		double *flux_faces;            // Riemann fluxes on the three interfaces
		double *quad_F1, *quad_F2;     // (F, G) at the quadrature points

		FE_2D_Basis* FE_BASIS;

		// Object that contains all the FE data needed
		FE_Store *FE_DATA;

		// Equation
		Equation *MODEL;

		// Element auxiliary methods
		virtual void compute_imass() = 0;
		virtual void Jacobian(Matrix2d &J, double xi, double eta) = 0;
		virtual void XI_to_X(double *x, double *y, double xi, double eta) = 0;

		void xy_at_node(double *x, double *y, unsigned int node);
                void save_dofs();
		void eval_at_quad(double* state, unsigned int n);
		
		// Test methods
		void print();  
	
		// Main methods
		void setup(int m, int p, double* points, FE_2D_Basis * fe_basis,
						FE_Store *fe_data, Equation *model);
		void compute_quad_fluxes();
		void assemble_rhs(double * RHS);
		void surface_integral(double * RHS);
		void volume_integral(double * RHS);
		void compute_CFL(double dt);
		void step(double dt);
		void step(double dt, double alpha1, double alpha2, double beta);

		Element();
		virtual ~Element();
	
};

class Tri_Element: public Element{

	public:
		Tri_Element(): Element(){NF = 3;};
	
		virtual void compute_imass();
		void Jacobian(Matrix2d &J, double xi, double eta);
		void XI_to_X(double* x, double *y, double xi, double eta);

};

class Quad_Element: public Element{

	public:
		Quad_Element(): Element(){NF = 4;};
	
		virtual void compute_imass();	
		virtual void Jacobian(Matrix2d &J, double xi, double eta);
		void XI_to_X(double* x, double *y, double xi, double eta);

};

#endif 
