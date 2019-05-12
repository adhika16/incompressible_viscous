#include "SIMPLE.hpp"

// =====================================================
// =====================================================
void SIMPLE::pressure_correction(std::vector<std::vector<double>> &p)
{
	// == reset class variables
	this->_p_n.resize(Pars::nx, std::vector<double>(Pars::ny, 0.0e0));
	this->_p_n_1.resize(Pars::nx, std::vector<double>(Pars::ny, 0.0e0));
	// == internal variables
	double a,b,c,d;
	double p_n_bar, dphi(0);
	int n_iter=0;
	// -- independent value of internal variables
	a = 2*(_dt/pow(_dx,2) + _dt/pow(_dy,2));
	b = -_dt/pow(_dx,2);
	c = -_dt/pow(_dy,2);

	// == calculating corrected pressure 
	do
	{
		for (int i = 1; i < _nx-1; i++) // excluding boundary
		{
			for (int j = 1; j < _ny-1; j++) // excluding boundary
			{
				// -- dependent value of internal variables
				d = Pars::rho/_dx*(_u_star[i][j] - _u_star[i-1][j]) + Pars::rho/_dy*(_v_star[i][j] - _v_star[i][j-1]);
				p_n_bar = -1/a*(b*_p_n[i+1][j] + b*_p_n[i-1][j] + c*_p_n[i][j+1] + c*_p_n[i][j-1] + d);
				_p_n_1[i][j] = _p_n[i][j] + Pars::omega*(p_n_bar - _p_n[i][j]);
				// printf("%f\n", _p_n_1[i][j]);
			}
		}
		n_iter += 1;
		get_error(dphi);
    cout <<  "\r" << dphi << " " << n_iter/*<< "\r"*/;
    // fflush ( stdin );
		// printf("%f\n", dphi);
		// -- update pressure correction
		_p_n = _p_n_1;
	} while(dphi > 1.0e-4);

}
// =====================================================
// =====================================================
void SIMPLE::get_error(double &dphi)
{
	// == reset error calculation
	dphi = 0.0e0;
	// == internal variables
	double sum_dp = 0.0e0;
	for (int i = 1; i < _nx-1; i++) // excluding boundary
	{
		for (int j = 1; j < _ny-1; j++) // excluding boundary
		{
			 dphi += std::abs(_p_n_1[i][j] - _p_n[i][j]);
			 sum_dp += std::abs(_p_n_1[i][j]);
		}
	}	
	dphi /= sum_dp;

}
// =====================================================
// =====================================================