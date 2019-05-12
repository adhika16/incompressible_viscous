#ifndef INCLUDED_GLOBAL
#include "global.hpp"
#endif

#ifndef INCLUDED_INITIALIZATION
#include "src/initialization.hpp"
#endif

#ifndef INCLUDED_SIMPLE
#include "src/SIMPLE.hpp"
#endif

int main(int argc, char const *argv[])
{
	// == dictionary
	std::vector<std::vector<double>> x_IJ; // x coordinate to store pressure
	std::vector<std::vector<double>> y_IJ; // y coordinate to store pressure
	std::vector<std::vector<double>> x_ij_u;  // x coordinate of forward staggered grid to store u
	std::vector<std::vector<double>> y_ij_u;  // y coordinate of forward staggered grid to store u
	std::vector<std::vector<double>> x_ij_v;  // x coordinate of forward staggered grid to store v
	std::vector<std::vector<double>> y_ij_v;  // y coordinate of forward staggered grid to store v
	std::vector<int> idx_in,idx_out; //Index Inlet and Outlet

	// Ghost Points
	std::vector<double> ul,ur;	//ul = points u left, ur = points u right
	std::vector<double> vt,vb;	//vt = points v top, vb = points v bottom
	// Properties
	std::vector<std::vector<double>> p;  // Pressure
	std::vector<std::vector<double>> u;  // Velocity in x-direction (u)
	std::vector<std::vector<double>> v;  // Velocity in y-direction (v)
	
	// === main ===
	printf("\n\t=================================================");
	printf("\n\t====== Incompressible Viscous Flows Solver ======");
	printf("\n\t=================================================\n");
	// == initialization ==
	Initialization _Initialization;
	_Initialization.init_domain(x_IJ, y_IJ, x_ij_u, y_ij_u, x_ij_v, y_ij_v,idx_in,idx_out,p,u,v);
	/*
	your code here ...
	*/
	// == SIMPLE method ==
	/*
	your code here ...
	*/
	SIMPLE _SIMPLE;
	_SIMPLE.get_SIMPLE(p, u, v, ul, ur, vt, vb, idx_in, idx_out, x_IJ, y_IJ, x_ij_u, y_ij_u, x_ij_v, y_ij_v);

	// == saving data ==
	// std::ofstream outs;
	// outs.open("output/airfoildata.dat");
	// for (int i = 0; i < x.size(); i++)
	// 	outs << x[i] << " " << y[i] << "\n";
	// outs.close();

	return EXIT_SUCCESS;
}
