#ifndef INCLUDED_GLOBAL
#include "global.hpp"
#endif

#ifndef INCLUDED_INITIALIZATION
#include "src/initialization.hpp"
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

	std::vector<std::vector<double>> p;  // Pressure
	std::vector<std::vector<double>> u;  // Velocity in x-direction (u)
	std::vector<std::vector<double>> v;  // Velocity in y-direction (v)
	// === main ===
	printf("\n\t=================================================");
	printf("\n\t====== Incompressible Viscous Flows Solver ======");
	printf("\n\t=================================================\n");
	// == initialization ==
	Initialization _Initialization;
	_Initialization.init_domain(x_IJ, y_IJ, x_ij_u, y_ij_u, x_ij_v, y_ij_v,p,u,v);
	/*
	your code here ...
	*/
	// == SIMPLE method ==
	/*
	your code here ...
	*/

	// == saving data ==
	// std::ofstream outs;
	// outs.open("output/airfoildata.dat");
	// for (int i = 0; i < x.size(); i++)
	// 	outs << x[i] << " " << y[i] << "\n";
	// outs.close();

	return EXIT_SUCCESS;
}
