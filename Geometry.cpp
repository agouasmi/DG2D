#include "Geometry.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

Geometry::Geometry(string mesh_name)
{
	ifstream file;

	string Positions = mesh_name +  ".Pos";
	string Elements  = mesh_name +  ".Ele";
	string Int_edges = mesh_name + ".IEdg";
	string Ext_edges = mesh_name + ".EEdg";

	// Read the positions
	file.open(Positions.c_str());
	file >> NODE_NUM;
	Pos = MatrixXd::Zero(NODE_NUM,2);

	int i = 0;
	double x,y;
	while (file >> x >> y){

		Pos(i,0) = x;
		Pos(i,1) = y;
		i++;	

	}
	file.close();
	
	// Read the elements
	file.open(Elements.c_str());
	file >> ELEM_NUM >> NF;
	Elem = MatrixXd::Zero(ELEM_NUM,NF);

	i = 0;
	int n1, n2, n3, n4;

	if (NF == 3){
		while (file >> n1 >> n2 >> n3){

			Elem(i,0) = n1-1; 
			Elem(i,1) = n2-1; 
			Elem(i,2) = n3-1;
			i++;	

		}
	}
	else {
       		while (file >> n1 >> n2 >> n3 >> n4){

			Elem(i,0) = n1-1; 
			Elem(i,1) = n2-1; 
			Elem(i,2) = n3-1;
			Elem(i,3) = n4-1;
			i++;	

		}
	}
	file.close();
		
	// Read the interior edges
	file.open(Int_edges.c_str());
	file >> I_EDGE_NUM;
	I_Edges = MatrixXd::Zero(I_EDGE_NUM,6);	

	i = 0;
	int v1, v2, E1, E2, F1, F2;
	int ele[3];

	while (file >> v1 >> v2 >> E1 >> E2){

		I_Edges(i,0) = v1-1; I_Edges(i,1) = v2-1;
		I_Edges(i,2) = E1-1; I_Edges(i,3) = E2-1;
		
		ele[0] = Elem(E1-1,0);
		ele[1] = Elem(E1-1,1);
		ele[2] = Elem(E1-1,2);

		F1 = track(ele, v1-1);
		
		ele[0] = Elem(E2-1,0);
		ele[1] = Elem(E2-1,1);
		ele[2] = Elem(E2-1,2);
		
		F2 = track(ele, v2-1);

		I_Edges(i,4) = F1; 
		I_Edges(i,5) = F2;
	
		i++;	

	}
	file.close();
	
	// Read the Exterior edges
	file.open(Ext_edges.c_str());
	file >> E_EDGE_NUM;
	E_Edges = MatrixXd::Zero(E_EDGE_NUM,4);

	i = 0;
	int E, F;
	while (file >> v1 >> v2 >> E){

		E_Edges(i,0) = v1-1; E_Edges(i,1) = v2-1;
		E_Edges(i,2) = E-1; 
		
		ele[0] = Elem(E-1,0);
		ele[1] = Elem(E-1,1);
		ele[2] = Elem(E-1,2);

		F = track(ele, v1-1);

		E_Edges(i,3) = F; 		

		i++;	

	}
	
	file.close();
}

int Geometry::track(int* tri, int node){

	int face;
		
	if (tri[0] == node){
		face = 0;
	}
	else{
		if (tri[1] == node){
			face = 1;
		}
		else{
			if (tri[2] == node){
				face = 2;
			}
			else{
				face = 3;
			}
		}
	}

	return face;

}

void Geometry::print(){

	cout << endl << "---- MESH ----" << endl << endl;
	cout << "There are " << NODE_NUM << " nodes" << endl;
	cout << "There are " << ELEM_NUM << " elements ( " << NF << " faces)" << endl;	
        cout << "There are " << I_EDGE_NUM << " Interior edges" << endl;
	cout << "There are " << E_EDGE_NUM << " Exterior edges" << endl;

}
/*
Triangulation::Triangulation(string mesh_name)
{
	ifstream file;

	string Positions = mesh_name +  ".Pos";
	string Triangles = mesh_name +  ".Tri";
	string Int_edges = mesh_name + ".IEdg";
	string Ext_edges = mesh_name + ".EEdg";

	// Read the positions
	file.open(Positions.c_str());
	file >> NODE_NUM;

	Pos = MatrixXd::Zero(NODE_NUM,2);

	int i = 0;
	double x,y;
	while (file >> x >> y){

		Pos(i,0) = x;
		Pos(i,1) = y;

		i++;	

	}
	
	file.close();
	
	// Read the Triangles
	file.open(Triangles.c_str());

	file >> TRI_NUM;

	Tri = MatrixXd(TRI_NUM,3);

	i = 0;
	int n1, n2, n3;

	while (file >> n1 >> n2 >> n3){

		Tri(i,0) = n1-1; Tri(i,1) = n2-1; Tri(i,2) = n3-1;
		i++;	

	}
	
	file.close();
		
	// Read the interior edges
	file.open(Int_edges.c_str());

	file >> I_EDGE_NUM;

	I_Edges = MatrixXd::Zero(I_EDGE_NUM,6);	

	i = 0;
	int v1, v2, T1, T2, F1, F2;
	while (file >> v1 >> v2 >> T1 >> T2 >> F1 >> F2){

		I_Edges(i,0) = v1-1; I_Edges(i,1) = v2-1;
		I_Edges(i,2) = T1-1; I_Edges(i,3) = T2-1;
		I_Edges(i,4) = F1; I_Edges(i,5) = F2;

		i++;	

	}
	
	file.close();
	
	// Read the Exterior edges
	file.open(Ext_edges.c_str());

	file >> E_EDGE_NUM;

	E_Edges = MatrixXd::Zero(E_EDGE_NUM,4);

	i = 0;
	int T, F;
	while (file >> v1 >> v2 >> T >> F){

		E_Edges(i,0) = v1-1; E_Edges(i,1) = v2-1;
		E_Edges(i,2) = T-1; 
		E_Edges(i,3) = F;

		i++;	

	}
	
	file.close();

	cout << endl;
}*/




