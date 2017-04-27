#include <iostream>
#include <stdio.h>
#include <string>
#include <fstream>

#include "DG_2D.h"

using namespace std;

/* --- Main program --- */

int main(int argc, char **argv){

	/* Read input parameters */
	
	double dt, T;
	int Nt, M, ORDER, freq, RK;
	string temp;

	ifstream file;
	file.open("Input/input.txt");
	file >> temp >> dt >> temp >> T >> temp >> M >> temp >> ORDER;
	file >> temp >> RK >> temp >> freq;
	file.close();

	Nt = int (T/dt);
	string mesh_path = "Mesh/Mesh";
	
	/* Problem Instantiation */
	
	Equation * MODEL = new Euler();

	DG_2D Problem(MODEL, ORDER, mesh_path);
	Problem.setup_system();

	Problem.print();
	Problem.write_sol("solution_0.txt");

	/* Main loop */

	cout << "---- RUN ----" << endl << endl;
	
	int nt = 0;
	int N_f = (int) Nt/freq;

	for (int f = 0; f < freq; f++){
		for (int i = 0; i < N_f; i++){
			Problem.advance_RK(dt,RK);
			nt++;	
		}
		Problem.write_sol("solution_" + to_string(f+1) + ".txt");
		
	}
	
	for (int i = nt; i < Nt; i++){
		Problem.advance_RK(dt,RK);
	}

	Problem.write_sol("solution_" + to_string(freq+1) + ".txt");
		
	delete MODEL;
	return 0;

}
