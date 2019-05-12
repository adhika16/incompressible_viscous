#ifndef INCLUDED_INITIALIZATION
#define INCLUDED_INITIALIZATION

#ifndef INCLUDED_GLOBAL
#include "../global.hpp"
#endif


class Initialization
{
  int _nx;
  int _ny;
	double _lx;
  double _ly;
  double _dx;
  double _dy;

public:
	Initialization();

  void init_domain(std::vector<std::vector<double>> &x_IJ, std::vector<std::vector<double>> &y_IJ,
   std::vector<std::vector<double>> &x_ij_u, std::vector<std::vector<double>> &y_ij_u,
   std::vector<std::vector<double>> &x_ij_v, std::vector<std::vector<double>> &y_ij_v,
   std::vector<int> &idx_in,std::vector<int> &idx_out,
   std::vector<std::vector<double>> &P,std::vector<std::vector<double>>&u, std::vector<std::vector<double>>&v);

	~Initialization();

};


#endif
