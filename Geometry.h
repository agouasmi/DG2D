#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <string>
#include </usr/include/eigen3/Eigen/Dense>

using namespace std;

using Eigen::MatrixXd;

/* Triangles and Quads */
class Geometry{

	public:
		int	NODE_NUM, NF, ELEM_NUM, I_EDGE_NUM, E_EDGE_NUM;  

		MatrixXd 		 Pos; 
		MatrixXd               	 Elem;

		MatrixXd    I_Edges, E_Edges; 

		Geometry(string mesh_file);
		int track(int* ele, int node);
		void print();		

};

#endif

