#ifndef INCLUDED_GLOBAL
#define INCLUDED_GLOBAL

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>
#include <algorithm>

using namespace std;

namespace Pars
{
	extern const double miu;
	extern const double rho;
	extern const double p_atm;
	extern const double lx;
	extern const double ly;
	extern const int nx;
	extern const int ny;
	extern const double dx;
	extern const double dy;
	extern const double omega;
	extern const double alpha_p;

	// Inlet and Outlet Properties
	extern const double p_i; //Inlet Pressure
	extern const double p_o; //Outlet Pressure
	extern const double v_i; //Inlet Velocity (v)
	extern const double u_i; //Inlet Velocity (u)

	// Inlet and Outlet Geometry
	extern const double l_i; //Inlet length
	extern const double l_o; //Outlet length

	// Time Step
	extern const double dt; 
}

#endif
