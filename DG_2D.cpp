#include "DG_2D.h"

#include <fstream>

using namespace std;

DG_2D::DG_2D(Equation* model, int p, string mesh_path){	
	
	M     = model->M;
	ORDER = p;

	MODEL = model;		
	MESH  = new Geometry(mesh_path);

	N_E = MESH->ELEM_NUM;
	t   = 0;

	QUAD1     = new Gauss(ORDER + 1);  
	BASIS     = new Lagrange(ORDER);
	
	if (MESH->NF == 3){ 
		ELEMENTS  = new Tri_Element[N_E];        
		QUAD2     = new Dunavant(2*ORDER + 1);       
		FE_BASIS  = new Lagrange_Tri(ORDER);
	}
	else{
		if (MESH->NF == 4){
			ELEMENTS = new Quad_Element[N_E];
			QUAD2    = new Tensor(QUAD1);
			FE_BASIS = new Lagrange_Quad(ORDER);
		}
	}
	
	FE_DATA = new FE_Store(QUAD1, QUAD2, FE_BASIS, BASIS);
	GLOBAL_DOF_NUM = N_E * FE_BASIS->DIM;

}

DG_2D::~DG_2D(){
	
	delete [] ELEMENTS;

	delete QUAD2;
	delete QUAD1;
	delete FE_BASIS;
	delete BASIS;
	
	delete FE_DATA;

	delete MESH;
}


void DG_2D::print(){

	//MESH->print();
	//QUAD1->print();
	//QUAD2->print();
	//FE_BASIS->print();
	//BASIS->print();
	//FE_DATA->print();

}

void DG_2D::setup_system(){

	FE_BASIS->setup();
	FE_DATA->compute(QUAD1, QUAD2, FE_BASIS, BASIS);
	
	// Initialize the ELEMENTS 
	int ELM;
	double points[2*MESH->NF];	
	
	for (ELM = 0; ELM < N_E; ELM++){

		for (int f = 0; f < MESH->NF; f++){
			points[2*f]   = MESH->Pos(MESH->Elem(ELM,f),0);
		        points[2*f+1] = MESH->Pos(MESH->Elem(ELM,f),1);
		}

		ELEMENTS[ELM].setup(M, ORDER, points, FE_BASIS, FE_DATA, MODEL);
	}

	set_initial_conditions();
	
}

void DG_2D::set_initial_conditions(){

	int ELEM, dof, m;
	double x, y, xi, eta, t0 = 0.;	

	for (ELEM = 0; ELEM < N_E; ELEM++){
		
		// Interpolation
		for (dof = 0; dof < FE_BASIS->DIM; dof++){
			ELEMENTS[ELEM].xy_at_node(&x, &y, dof);
			MODEL->init_sol(&(ELEMENTS[ELEM].dofs[dof*M]), x, y);
		}

		// Least-squares projection
		/*MatrixXd B(FE_BASIS->DIM,M);
		Map < Matrix<double, Dynamic, Dynamic, RowMajor> > 
					U(ELEMENTS[ELEM].dofs, FE_BASIS->DIM, M);

		double state[M], temp[M], J, w;
		J = ELEMENTS[ELEM].J;
	
		for (dof = 0; dof < FE_BASIS->DIM; dof++){

			for (m = 0; m < M; m++){
				state[m] = 0;
			}

			for (int n = 0; n < QUAD2->N_PTS; n++){

				xi  = QUAD2->POINTS[2*n];
				eta = QUAD2->POINTS[2*n+1];
				
				w = QUAD2->WEIGHTS[n];
	
				ELEMENTS[ELEM].XI_to_X(&x, &y, xi, eta);
				MODEL->init_sol(temp, x, y);
			
				for (m = 0; m < M; m++){
					state[m] += w*temp[m]*FE_BASIS->PHI(dof,xi,eta);
				}

			}
			
			for ( m = 0; m < M; m++){
				B(dof,m) = J*state[m];
			}
		}
		
		U = ELEMENTS[ELEM].IMASS*B;*/
	
	}

}


void DG_2D::get_CFL(){
	
	/*double temp = ELEMENTS[0].CFL;
	for (int ELEM = 1; ELEM < N_E; ELEM++){
		
		if (ELEMENTS[ELEM].CFL > temp){
			temp = ELEMENTS[ELEM].CFL;
		}

	}
	CFL = temp;*/
}

void DG_2D::apply_BC(){	
}

void DG_2D::compute_fluxes(){
	/* Compute the fluxes at the interior edges */

	int EDGE, Ele, Ele_L, Ele_R, Face, Face_L, Face_R, point;
	int ept, ept_L, ept_R, dof, left_dof, right_dof;
	double x, y, nx, ny;	
	double smax;

	for (EDGE = 0; EDGE < MESH->I_EDGE_NUM; EDGE++){

		//smax = 0;

		Ele_L  = MESH->I_Edges(EDGE, 2);
		Ele_R  = MESH->I_Edges(EDGE, 3);
	
		Face_L = MESH->I_Edges(EDGE, 4);
		Face_R = MESH->I_Edges(EDGE, 5); 

		nx  = ELEMENTS[Ele_L].N(Face_L,0);
		ny  = ELEMENTS[Ele_L].N(Face_L,1);
				
		for (point = 0; point < ORDER + 1; point++){
			
			ept_L = Face_L*(ORDER + 1) + point;
			ept_R = Face_R*(ORDER + 1) + ORDER - point;

			left_dof  = FE_BASIS->face_index[ept_L];
			right_dof = FE_BASIS->face_index[ept_R];

			//smax += 
			MODEL->interface_flux(&(ELEMENTS[Ele_L].flux_faces[ept_L*M]), 
				         &(ELEMENTS[Ele_R].flux_faces[ept_R*M]),
				         nx,ny, 
				         &(ELEMENTS[Ele_L].dofs[left_dof*M]),
				         &(ELEMENTS[Ele_R].dofs[right_dof*M]));
		}

		//smax /= ORDER + 1;
		//ELEMENTS[Tri_L].smax[Face_L] = smax;
		//ELEMENTS[Tri_R].smax[Face_R] = smax;
	}
        
	for (EDGE = 0; EDGE < MESH->E_EDGE_NUM; EDGE++){

		//smax = 0;

		Ele  = MESH->E_Edges(EDGE, 2);
		Face = MESH->E_Edges(EDGE, 3);
		
		nx  = ELEMENTS[Ele].N(Face,0);
		ny  = ELEMENTS[Ele].N(Face,1);
				
		for (point = 0; point < ORDER + 1; point++){
			
			ept = Face*(ORDER + 1) + point;
			dof = FE_BASIS->face_index[ept];

			ELEMENTS[Ele].xy_at_node(&x, &y, dof);		

			//MODEL->BC_flux(&(ELEMENTS[Ele].flux_faces[ept*M]),
				// x, y, nx, ny, &(ELEMENTS[Ele].dofs[dof*M]));
		}

	}
	//get_CFL();
}

void DG_2D::advance(double dt){
		
	int ELEM;
	
	save_dofs();

	apply_BC();
	compute_fluxes();

	for (ELEM = 0; ELEM < N_E; ELEM++){
		ELEMENTS[ELEM].step(dt);
	}

	t += dt;
}

void DG_2D::save_dofs(){

	for (int ELEM = 0; ELEM < N_E; ELEM++){
		ELEMENTS[ELEM].save_dofs();
	}

}

void DG_2D::advance_RK(double dt, int RK){

	int ELEM, level;

	save_dofs();
	
	for (level = RK; level > 0; level--){
		
		apply_BC();
		compute_fluxes();
	
		for (ELEM = 0; ELEM < N_E; ELEM++){
			ELEMENTS[ELEM].step(dt/(double) level);
		}
	}

	t += dt;
}

void DG_2D::print_elements(){

	int ELEM;
	for (ELEM = 0; ELEM < N_E; ELEM++){

		cout << " Element " << ELEM << ": " << endl;
		ELEMENTS[ELEM].print();
	} 
}

void DG_2D::write_sol(string filename){
	
	string out = "Output/" + filename;	

	ofstream myfile;
  	myfile.open (out.c_str());

  	int ELEM, node, m, r, s;
	double x, y;
	
	for (ELEM = 0; ELEM < N_E; ELEM++){
		for (node = 0; node < FE_BASIS->DIM; node++){
		      	
			ELEMENTS[ELEM].xy_at_node(&x, &y, node);
			myfile << x << " " << y ;

			for (int m = 0; m < M; m++){
					myfile << " " << ELEMENTS[ELEM].dofs[node*M+m];
			}		

			myfile << endl;
		}	
	}
	
        myfile.close();
}
