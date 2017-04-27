/*       -------- 2D Discontinuous Galerkin Class  --------- */

#ifndef DG_2D_H
#define DG_2D_H

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string>

#include "Element.h"
#include "Equation.h"
#include "FE_Store.h"
#include "Geometry.h"

using namespace std;

class DG_2D{
	
	public:
		DG_2D(Equation* model, int p, string mesh_path);
                ~DG_2D();
	
		double t;
		double CFL;

		Equation* MODEL;
		Geometry*  MESH;           // Has all the information you want
	
		int    M, ORDER, N_E;          // m: number of equations
		int    GLOBAL_DOF_NUM;         // total number of dofs = (N_E * DIM)

		Quadrature_2D * QUAD2;
		Quadrature_1D * QUAD1;            
		FE_2D_Basis   * FE_BASIS;
		Basis_1D      * BASIS;         

		FE_Store * FE_DATA;

		Element * ELEMENTS;  
	
		// Solution steps
		void apply_BC();
		void compute_fluxes();
		void get_CFL();
		void advance(double dt);
		void advance_RK(double dt, int RK);
		
		// Auxiliary functions
		void print();
                void print_elements();
		void write_sol(string);		
		void save_dofs();	
	
		void setup_system();           
		void set_initial_conditions();		
};


#endif
