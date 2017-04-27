#include "Element.h"

/* Base class */
Element::Element(){}

Element::~Element(){

	//delete [] smax;

	// ---------------- //

	delete [] dofs_rhs;
	delete [] dofs;
	delete [] dofs_backup;
	
	// ---------------- //	
	
	delete [] quad_F1; 
	delete [] quad_F2;

	// ---------------- //
	
	delete [] flux_faces; 
	
	// --------------- //

	FE_BASIS = NULL;
	FE_DATA  = NULL;
	MODEL    = NULL;	
}

void Element::setup(int m, int p, double* points, FE_2D_Basis* fe_basis,
					FE_Store *fe_data, Equation* model){
	
	// Parameters
	M           = m;
	ORDER       = p;
	DOFS_NUMBER = fe_data->FE_DIM;

	FE_BASIS = fe_basis;
	MODEL = model;
	FE_DATA = fe_data;

	N = MatrixXd(NF,2);
	X = VectorXd(NF,1);
	Y = VectorXd(NF,1);
	L = VectorXd(NF,1);

	// Geometry info
	for (int i = 0; i < NF; i++){
		X(i) = points[2*i];
	        Y(i) = points[2*i+1];
	} 
	
	for (int i = 0; i < NF-1; i++){
		N(i,0) = Y(i+1) - Y(i); 
		N(i,1) = - (X(i+1) - X(i));		
	}
	N(NF-1,0) =    Y(0) - Y(NF-1);
	N(NF-1,1) = - (X(0) - X(NF-1));
	
	//smax = new double[NF];
	double P = 0.;
	for (int i = 0; i < NF; i++){
		L(i) = N.row(i).norm();
		N.row(i) /= L(i);
		P += L(i); 
	}

	size = 1.0;
	
	compute_imass();	
	
	// Allocating memory
	dofs_rhs    = new double[M*DOFS_NUMBER];
	dofs        = new double[M*DOFS_NUMBER];
	dofs_backup = new double[M*DOFS_NUMBER];

	quad_F1 = new double[M*FE_DATA->Q2_N_PTS];
	quad_F2 = new double[M*FE_DATA->Q2_N_PTS];

	flux_faces = new double[NF*M*(ORDER+1)]; 
	
}

void Element::xy_at_node(double *x, double *y, unsigned int node){
	
	double xi, eta;

	xi  = FE_DATA->XIETA[2*node];
	eta = FE_DATA->XIETA[2*node + 1];

	XI_to_X(x, y, xi, eta);
}	


void Element::eval_at_quad(double * U, unsigned int n){
	/* Evaluates the state in the reference element */

	double phi;

	for (int m = 0; m < M; m++){
		U[m] = 0.;
	}
		
	for (int dof = 0; dof < DOFS_NUMBER; dof++){
		phi = FE_DATA->PHI[dof*(FE_DATA->Q2_N_PTS) + n];
		for (int m = 0; m < M; m++){
			U[m] += dofs[dof*M + m]*phi;
		}
	}

}


void Element::compute_quad_fluxes(){
	/* Computes the weighted flux values at 2D quadrature points */	
	
	double quad_dof[M], xi, eta, w;
	int n, m;
	
	for (n = 0; n < FE_DATA->Q2_N_PTS; n++){

		w = FE_DATA->Q2_WEIGHTS[n];
		
		// Eval the state at the quad point
		eval_at_quad(quad_dof, n);

		MODEL->flux_N(&quad_F1[n*M], quad_dof, 1., 0.);
		MODEL->flux_N(&quad_F2[n*M], quad_dof, 0., 1.);
	
		// Multiply by the weight	
		for (m = 0; m < M; m++){
			quad_F1[n*M + m] *= w;
			quad_F2[n*M + m] *= w;
		}

	}
}

void Element::volume_integral(double* RHS){

	int node, n;
	double xi, eta, dphi_dxi, dphi_deta, detJ;

	RowVector2d NABLA_PHI;
	MatrixXd F(2,M);
	Matrix2d J_, J1;

	Map < RowVectorXd > F1(NULL,M);
	Map < RowVectorXd > F2(NULL,M);
	Map < RowVectorXd > R_node(NULL,M);
	
	for (node = 0; node < DOFS_NUMBER; node++){

		new (&R_node) Map < RowVectorXd >(&RHS[node*M],M);		
		R_node *= 0.;
		
		for (n = 0; n < FE_DATA->Q2_N_PTS; n++){
	
			new (&F1) Map < RowVectorXd > (&quad_F1[n*M],M);
			new (&F2) Map < RowVectorXd > (&quad_F2[n*M],M);
			
			F.row(0) = F1;
			F.row(1) = F2;

			xi  = FE_DATA->Q2_POINTS[2*n];
			eta = FE_DATA->Q2_POINTS[2*n+1];
	
			Jacobian(J_, xi, eta); 
			detJ = J_.determinant();
			J1 = J_.inverse();

			dphi_dxi  = FE_DATA->DPHI_DXI [node*FE_DATA->Q2_N_PTS + n];
			dphi_deta = FE_DATA->DPHI_DETA[node*FE_DATA->Q2_N_PTS + n];

			NABLA_PHI << dphi_dxi , dphi_deta;

			R_node += detJ*NABLA_PHI*J1*F;
		}
		
	}	

}

void Element::surface_integral(double* RHS){

	Map < RowVectorXd > F_edge(NULL,M);
	Map < RowVectorXd > R_node(NULL,M);	
	
	int face, k, node, n;
	double w;
	
	for (face = 0; face < NF; face++){
		for (k = 0; k < ORDER + 1; k++){

			node = FE_DATA->face_index[(ORDER+1)*face + k];

			new (&R_node) Map < RowVectorXd > (&RHS[node*M],M);
			new (&F_edge) Map < RowVectorXd > (&flux_faces[((ORDER+1)*face + k)*M],M);
		
			for (n = 0; n < FE_DATA->Q1_N_PTS; n++){
							
				w  = FE_DATA->Q1_WEIGHTS[n];
				R_node -= 0.5 * L[face] * w * F_edge *
					   FE_DATA->EDGE_PHI[k*(FE_DATA->Q1_N_PTS) + n];			
			}
	
		} 
	}
}

void Element::assemble_rhs(double * RHS){

	// Volume term
	volume_integral(RHS);
	
	// Interface term
	surface_integral(RHS);		

}

void Element::compute_CFL(double dt){

	/*double S_MAX, P = J/size;
	for (int face = 0; face < NF; face++){
		S_MAX += smax[face]*L[face]/P;
	}

	CFL = dt * S_MAX / size;*/
}

void Element::step(double dt){
	/* Euler explicit */
	
	compute_quad_fluxes();	

	assemble_rhs(dofs_rhs);
	
	Map< Matrix<double, Dynamic, Dynamic, RowMajor> > RHS(dofs_rhs, DOFS_NUMBER, M);
	Map< Matrix<double, Dynamic, Dynamic, RowMajor> > U(dofs, DOFS_NUMBER, M);
	Map< Matrix<double, Dynamic, Dynamic, RowMajor> > U_OLD(dofs_backup, DOFS_NUMBER, M);
	compute_CFL(dt);

	U = U_OLD + dt*IMASS*RHS;	
}

void Element::step(double dt, double alpha1, double alpha2, double beta){
	/* SSP Euler explicit */
	
	compute_quad_fluxes();	

	assemble_rhs(dofs_rhs);
	
	Map< Matrix<double, Dynamic, Dynamic, RowMajor> > RHS(dofs_rhs, DOFS_NUMBER, M);
	Map< Matrix<double, Dynamic, Dynamic, RowMajor> > U(dofs, DOFS_NUMBER, M);
	Map< Matrix<double, Dynamic, Dynamic, RowMajor> > U_OLD(dofs_backup, DOFS_NUMBER, M);

	U = alpha1*U_OLD + alpha2*U + dt*IMASS*RHS*beta;	
}

void copy_vec(double *out, double* in, int N){
	memcpy(out, in, N * sizeof (double));
}

void Element::save_dofs(){
	copy_vec(dofs_backup, dofs, M*DOFS_NUMBER);
}

void Element::print(){
	/* Prints Lagrangian-like information on the element */

	double x, y;

	cout << " Vertices positions: " << X << endl << endl;
	cout << " Positions of the nodes:" << endl;

	for (int node = 0; node < DOFS_NUMBER; node++){		
		
		xy_at_node(&x, &y, node);
		cout << " Node " << node << ": ";
		cout << x << "  " << y << endl; 
	}

	cout << endl << " ----------------- " << endl;
}

/* Triangle */
void Tri_Element::compute_imass(){
	
	Matrix2d J;
	Jacobian(J, 0., 0.);
	 
	IMASS = FE_DATA->MASS_INV / J.determinant();
	
}

void Tri_Element::Jacobian(Matrix2d &J, double xi, double eta){
	
	J << X(1) - X(0), X(2) - X(0),
	     Y(1) - Y(0), Y(2) - Y(0);
		
}

void Tri_Element::XI_to_X(double *x, double *y, double xi, double eta){

	*x = X(0) + (X(1) - X(0))*xi + (X(2) - X(0))*eta;
	*y = Y(0) + (Y(1) - Y(0))*xi + (Y(2) - Y(0))*eta;

}

/* Quadrilateral */
void Quad_Element::compute_imass(){

	Matrix2d J_;
	MatrixXd MASS(DOFS_NUMBER, DOFS_NUMBER);

	int i, j, n;
	double xi, eta, detJ;

	Gauss quad1(ORDER+1);
	Tensor Quad(&quad1);
	
	for (i = 0; i < DOFS_NUMBER; i++){
		for (j = 0; j < DOFS_NUMBER; j++){
			
			MASS(i,j) = 0.;
			for (n = 0; n < Quad.N_PTS; n++){

				xi  = Quad.POINTS[2*n];
				eta = Quad.POINTS[2*n+1];

				Jacobian(J_, xi, eta);
				detJ = J_.determinant();

				MASS(i,j) += detJ*Quad.WEIGHTS[n]*
						FE_BASIS->PHI(i,xi,eta)*
						FE_BASIS->PHI(j,xi,eta);

			}			

		}	
	}
	
	IMASS = MASS.inverse();
}

void Quad_Element::Jacobian(Matrix2d &J, double xi, double eta){
		
	J(0,0) = 0.25*((X(1) - X(0))*(1 - eta) + (X(2) - X(3))*(1 + eta));
	J(1,0) = 0.25*((Y(1) - Y(0))*(1 - eta) + (Y(2) - Y(3))*(1 + eta));
	J(0,1) = 0.25*((X(3) - X(0))*(1 - xi)  + (X(2) - X(1))*(1 + xi));
	J(1,1) = 0.25*((Y(3) - Y(0))*(1 - xi)  + (Y(2) - Y(1))*(1 + xi));

}

void Quad_Element::XI_to_X(double *x, double *y, double xi, double eta){

	*x = 0.25*(X(0)*(1-xi)*(1-eta) + X(1)*(1+xi)*(1-eta) 
			+ X(2)*(1+xi)*(1+eta) + X(3)*(1-xi)*(1+eta));
	*y = 0.25*(Y(0)*(1-xi)*(1-eta) + Y(1)*(1+xi)*(1-eta) 
			+ Y(2)*(1+xi)*(1+eta) + Y(3)*(1-xi)*(1+eta));

}
