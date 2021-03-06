#include "global.hpp"

namespace Pars
{
	extern const double miu = 1.0e-5;
	extern const double rho = 1.226;
	extern const double p_atm = 1.0e5;
	extern const double lx = 1;
	extern const double ly = 1;
	extern const int nx = 50;
	extern const int ny = 50;
	extern const double dx = lx/(double)nx;
	extern const double dy = ly/(double)ny;
	extern const double omega = 1.5; // relaxation factor
	extern const double alpha_p = 0.01; // pressure correctoin relaxation factor

	// Inlet and Outlet Properties
	extern const double p_i = 0.0; //Inlet Pressure
	extern const double p_o = 0.0; //Outlet Pressure
	extern const double v_i = 0.0; //Inlet Velocity (v)
	extern const double u_i = 1.0; //Inlet Velocity (u)
	
	// Inlet and Outlet Geometry
	extern const double l_i = 0.25; //Inlet length
	extern const double l_o = 0.25; //Outlet length

	// Time Step
	extern const double dt = 0.0001;
}
