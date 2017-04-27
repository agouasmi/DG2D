#include "Equation.h"
#include <stdio.h>
#include <string.h>

using namespace std;

/* Euler class */

void prim2conv(double* conv, double* prim);
void conv2prim(double* prim, double* conv);
void exact_sol(double* conv, double x, double y, double t);

void prim2conv(double *conv, double* prim){
	
	/* output: r, r*u, r*v, r*E */

	double r, u, v, p;
	
	r = prim[0]; 
	u = prim[1]; 
	v = prim[2]; 
	p = prim[3];

	conv[0] = r;
	conv[1] = r*u;
	conv[2] = r*v;
	conv[3] = p/(GAMMA-1) + 0.5*r*(pow(u,2) + pow(v,2));
	
}

void conv2prim(double *prim, double* conv){

	/* output: rho, u, v, p, H */
	
	double r, ru, rv, rE;
	
	r  = conv[0];
	ru = conv[1];
	rv = conv[2];
	rE = conv[3];

	double p, H;

	p = (GAMMA - 1)*(rE - 0.5/r*(pow(ru,2) + pow(rv,2)));
	H = (rE + p)/r;

	prim[0] = r;
	prim[1] = ru/r;
	prim[2] = rv/r;
	prim[3] = p;
	prim[4] = H; 

}

void exact_sol(double* conv, double x, double y, double t){

	double u_inf, v_inf, prim[4], x0, y0, Beta, r, du, dv, T;

	x0 = 0.; y0 = 0.;
	u_inf = 1./sqrt(2.); v_inf = 1./sqrt(2.);
	Beta = 5;

	r = sqrt(pow(x - x0, 2) + pow(y - y0, 2));

	du = Beta/(2*M_PI)*exp(0.5*(1 - r*r))*(-(y - y0));
	dv = Beta/(2*M_PI)*exp(0.5*(1 - r*r))*( (x - x0));

	T = 1 - (GAMMA-1)*Beta*Beta/(8*GAMMA*M_PI*M_PI)*exp(1-r*r);

	prim[0] = pow(T,1./(GAMMA - 1.));
	prim[1] = u_inf + du;
	prim[2] = v_inf + dv;
	prim[3] = prim[0]*T;

	prim2conv(conv, prim);

}

void Euler::init_sol(double* conv, double x, double y){
	exact_sol(conv, x, y, 0.);
}

void Euler::flux_N(double * flux, double * conv, double nx, double ny){

	double prim[5], r, u, v, p, H;
	conv2prim(prim,conv);

	r = prim[0]; u = prim[1]; v = prim[2]; p = prim[3]; H = prim[4];

	flux[0] = nx * (r*u)            + ny * (r*v);
	flux[1] = nx * (r*pow(u,2) + p) + ny * (r*u*v);
	flux[2] = nx * (r*u*v)          + ny * (r*pow(v,2) + p);
	flux[3] = nx * (r*u*H)          + ny * (r*v*H);

}

void Euler::interface_flux(double *out1, double *out2, double nx, double ny,
						double * UL, double * UR){
	Roe_flux(out1, out2, nx, ny, UL, UR);
}

double Euler::Roe_flux(double * F, double * F_, double nx, double ny, 
				double * UL, double * UR)
{	double GMI = GAMMA - 1.0;

	// Process left state
	double rL, uL, vL, unL, qL, pL, rHL, HL, cL;

	rL  = UL[0];
	uL  = UL[1]/rL;
	vL  = UL[2]/rL;
	unL = uL*nx + vL*ny;
	qL = sqrt(pow(UL[1],2) + pow(UL[2],2))/rL;
	pL = (GMI)*(UL[3] - 0.5*rL*pow(qL,2));
	
	if ( (pL <= 0) || (rL <=0))
		cout << "Non-physical state!" << endl;
	
	rHL = UL[3] + pL;
	HL = rHL/rL;
	cL = sqrt(GAMMA*pL/rL);

	// Process right state
	double rR, uR, vR, unR, qR, pR, rHR, HR, cR;

	rR  = UR[0];
	uR  = UR[1]/rR;
	vR  = UR[2]/rR;
	unR = uR*nx + vR*ny;
	qR = sqrt(pow(UR[1],2) + pow(UR[2],2))/rR;
	pR = (GMI)*(UR[3] - 0.5*rR*pow(qR,2));
	
	if ( (pR <= 0) || (rR <=0))
		cout << "Non-physical state !" << endl;
	
	rHR = UR[3] + pR;
	HR = rHR/rR;
	cR = sqrt(GAMMA*pR/rR);

	// Left and right fluxes
	double FL[4] = {rL*unL, UL[1]*unL + pL*nx, UL[2]*unL + pL*ny, rHL*unL};
	double FR[4] = {rR*unR, UR[1]*unR + pR*nx, UR[2]*unR + pR*ny, rHR*unR};

	// Difference in states
	double du[4] = {UR[0] - UL[0], UR[1] - UL[1], UR[2] - UL[2], UR[3] - UL[3]};

	// Roe average
	double di, d1, ui, vi, Hi, af, ucp, c2, ci, ci1;
	
	di = sqrt(rR/rL);
	d1 = 1.0/(1.0+di);
	
	ui = (di*uR + uL)*d1;
	vi = (di*vR + vL)*d1;
	Hi = (di*HR + HL)*d1;

	af  = 0.5*(ui*ui + vi*vi);
	ucp = ui*nx + vi*ny;
	c2  = GMI*(Hi - af);
	
	if (c2<=0)
		cout << "Non-physical state !" << endl;

	ci  = sqrt(c2);
	ci1 = 1.0/ci;
	
	// Eigenvalues
	double l[3] = {ucp + ci, ucp - ci, ucp};
	
	// Entropy fix
	double epsilon = ci*.1;
	for (int i = 0; i < 3; i++){
		if ( abs(l[i]) < epsilon ) {
			l[i] = 0.5*(epsilon + l[i]*l[i]/epsilon);
		}
		l[i] = abs(l[i]);
	}	

	double l3 = l[2];

	double s1, s2, G1, G2, C1, C2;

	s1 = 0.5*(l[0] + l[1]);
	s2 = 0.5*(l[0] - l[1]);

	G1 = GMI*(af*du[0] - ui*du[1] - vi*du[2] + du[3]);
	G2 = -ucp*du[0] + du[1]*nx + du[2]*ny;
	
	C1 = G1*(s1 - l3)*ci1*ci1 + G2*s2*ci1;
	C2 = G1*s2*ci1            + G2*(s1 - l3);

	double smag = l[0];
	if (smag < l[1])
		smag = l[1];
	if (smag < l[2])
		smag = l[2];

	F[0]    = 0.5*(FL[0] + FR[0]) - 0.5*(l3*du[0] + C1);
	F[1]    = 0.5*(FL[1] + FR[1]) - 0.5*(l3*du[1] + C1*ui + C2*nx);
	F[2]    = 0.5*(FL[2] + FR[2]) - 0.5*(l3*du[2] + C1*vi + C2*ny);
	F[3]    = 0.5*(FL[3] + FR[3]) - 0.5*(l3*du[3] + C1*Hi + C2*ucp);
        
	/*F[0]    = 0.5*(FL[0] + FR[0]) - 0.5*smag*du[0];
	F[1]    = 0.5*(FL[1] + FR[1]) - 0.5*smag*du[1];
	F[2]    = 0.5*(FL[2] + FR[2]) - 0.5*smag*du[2];
	F[3]    = 0.5*(FL[3] + FR[3]) - 0.5*smag*du[3];*/

	for (int i = 0; i < M; i++){
		F_[i] = -F[i];
	}
	
	return smag;
}

void Euler::BC_flux(double * out, double x, double y, double nx, double ny, double *U){

	//double Ud[M] = {0., 0., 0., 0.};
	//flux_N(out, U, nx, ny);
}

/* Advection class */ 

void Advection::init_sol(double* conv, double x, double y){

	double t = 0.;
        double x0 = 0., y0 = 0.;
        double V_inf[2] = {VX, VY}, V_INF = sqrt(VX*VX + VY*VY);
        double rc    = 1.0, eps  = 0.3, M_inf = 0.5, r_inf = 1.;

        double e_f0, f1, f2;

        e_f0 = exp(1. - 1./pow(rc,2)*( pow(x - x0 - V_inf[0]*t, 2) +
                                                 pow(y - y0 - V_inf[1]*t, 2) ) );
        f1   = 1 - pow(eps,2)*(GAMMA - 1)*pow(M_inf,2)/(8.*pow(M_PI,2))*e_f0;
        f2   = eps * V_INF/(2*M_PI*rc)*(sqrt(e_f0));

        conv[0] = r_inf * pow(f1, 1./(GAMMA - 1));
        conv[1] = V_inf[0] - f2*(y - y0 - V_inf[1]*t);
}

void Advection::flux_N(double * flux, double * conv, double nx, double ny){

        double r = conv[0], u = conv[1];

        flux[0] = nx * (VX * r) + ny * (VY * r);
        flux[1] = nx * (VX * u) + ny * (VY * u);

}

void Advection::interface_flux(double *out1, double *out2, double nx, double ny,
						double * UL, double * UR){
	Roe_flux(out1, out2, nx, ny, UL, UR);
}

double Advection::Roe_flux(double * F, double * F_, double nx, double ny,
                                double * UL, double * UR)
{       // Computes the Roe flux (basic upwinding here)

        F[0] = 0.5 * (VX * nx + VY * ny) * (UL[0] + UR[0])
                        - 0.5 * abs(VX * nx + VY * ny) * (UR[0] - UL[0]);
        F[1] = 0.5 * (VX * nx + VY * ny) * (UL[1] + UR[1])
                        - 0.5 * abs(VX * nx + VY * ny) * (UR[1] - UL[1]);

        for (int i = 0; i < 2; i++){
                F_[i] = -F[i];
        }
}

void Advection::BC_flux(double * out, double x, double y, double nx, double ny, double *U){

	double Ud[M] = {0., 0., 0., 0.};
	flux_N(out, Ud, nx, ny);
}
