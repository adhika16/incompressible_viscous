#ifndef INCLUDED_SIMPLE
#define INCLUDED_SIMPLE

#ifndef INCLUDED_GLOBAL
#include "../global.hpp"
#endif


class SIMPLE
{
	int _nx,_ny;
	double _dx,_dy,_dt;
	std::vector<std::vector<double>> _u_star;
	std::vector<std::vector<double>> _v_star;

	std::vector<std::vector<double>> _p_n; // p'^{n}
	std::vector<std::vector<double>> _p_n_1;// p'^{n+1}

	void get_error(double &dphi);

public:
	SIMPLE();

  void MomentumEq(const std::vector<std::vector<double>> p,const std::vector<std::vector<double>> u,
    const std::vector<std::vector<double>> v,std::vector<double> ul,std::vector<double> ur,
    std::vector<double> vt,std::vector<double> vb,const std::vector<int> idx_in,
    const std::vector<int> idx_out);

	void pressure_correction(std::vector<std::vector<double>> &p);

	void get_SIMPLE(std::vector<std::vector<double>> p, std::vector<std::vector<double>> u,
    std::vector<std::vector<double>> v,std::vector<double> ul,std::vector<double> ur,
    std::vector<double> vt,std::vector<double> vb,const std::vector<int> idx_in,
    const std::vector<int> idx_out, const std::vector<std::vector<double>> x_IJ,
		const std::vector<std::vector<double>> y_IJ,
	  const std::vector<std::vector<double>> x_ij_u, const std::vector<std::vector<double>> y_ij_u,
	  const std::vector<std::vector<double>> x_ij_v, const std::vector<std::vector<double>> y_ij_v);

};


#endif
