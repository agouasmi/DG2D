#ifndef EQUATION_H
#define EQUATION_H

#include <math.h>
#include <cmath>
#include <iostream>

class Equation {

	public:
		int M; // number of states

		virtual void init_sol(double *conv, double x, double y) = 0;
		virtual void flux_N(double *flux, double *conv, 
						double nx, double ny) = 0;
		virtual void interface_flux(double *out1, double * out2,
						double nx, double ny,
						double* UL, double* UR) = 0;
		virtual void BC_flux(double * flux, double x, double y,
						double nx, double ny, double* U) = 0;

		Equation(int m):M(m) {};
		virtual ~Equation(){};

};

class Euler: public Equation {

	public:
	
		void init_sol(double* conv, double x, double y);
		void flux_N(double* flux, double *conv, double nx, double ny);
		void interface_flux(double* out1, double* out2, double nx, double ny,
							double *UL, double *UR);
		void BC_flux(double * flux, double x, double y, double nx, double ny,
									double * U);

		double Roe_flux(double* out1, double* out2, double nx, double ny,
								double* UL, double* UR);	
		Euler(): Equation(4){};		

};

class Advection: public Equation {

	public:
		double VX, VY;
		
		void init_sol(double* conv, double x, double y);
		void flux_N(double* flux, double *conv, double nx, double ny);
		void interface_flux(double* out1, double* out2, double nx, double ny,
							double *UL, double *UR);
		void BC_flux(double * flux, double x, double y, double nx, double ny,
									double * U);

		double Roe_flux(double* out1, double* out2, double nx, double ny,
								double* UL, double* UR);	
		Advection(double vx, double vy): Equation(2), VX(vx), VY(vy) {};		

};

#define GAMMA 1.4

#endif
