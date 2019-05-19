#include "SIMPLE.hpp"

void SIMPLE::get_SIMPLE(std::vector<std::vector<double>> &p, std::vector<std::vector<double>> &u,
  std::vector<std::vector<double>> &v, std::vector<double> ul,std::vector<double> ur,
  std::vector<double> vt,std::vector<double> vb,const std::vector<int> idx_in, const std::vector<int> idx_out, 
  const std::vector<std::vector<double>> &x_IJ, const std::vector<std::vector<double>> &y_IJ,
  const std::vector<std::vector<double>> &x_ij_u, const std::vector<std::vector<double>> &y_ij_u,
  const std::vector<std::vector<double>> &x_ij_v, const std::vector<std::vector<double>> &y_ij_v
	)
{
	// == internal variables
	double err_p, err_u, err_v;
	int n_iter = 0;
	printf("\n");
	do
	{
		this->MomentumEq(p, u, v, ul, ur, vt, vb, idx_in, idx_out);
		this->pressure_correction(p);
	// 	// == update primitive variables using converged pressure correction
		err_p=0;
		err_u=0;
		err_v=0;
		
		for (int i = 1; i < _nx-1; i++) // excluding boundary
		{
			for (int j = 1; j < _ny-1; j++) // excluding boundary
			{
				// p[i][j] += Pars::alpha_p*_p_n[i][j];
				// err_p   += abs(Pars::alpha_p*_p_n[i][j]);
				p[i][j] = p[i][j] + Pars::alpha_p*_p_n(i,j);
				err_p   += abs(Pars::alpha_p*_dx*_dy*_p_n(i,j));
			}
		}
		n_iter+=1;
		// printf("\n%f\n", err_p);
		// printf("\n");
    // cout <<  "\rrelative error:" << err_p << ", no. iterations: " << n_iter/*<< "\r"*/;

		u = this->_u_star;
		v = this->_v_star;

		for (int i = 0; i < _nx-2; i++) // excluding boundary
		{
			for (int j = 0; j < _ny; j++) // excluding boundary
			{
				// u[i][j] += -_dt/_dx*(_p_n[i+1][j] - _p_n[i][j]);
				// err_u   += abs(_dt/_dx*(_p_n[i+1][j] - _p_n[i][j]));
				// u[i][j] = _u_star[i][j] - _dt/(_dx*Pars::rho)*(_p_n(i+1,j) - _p_n(i,j));
				// u[i][j] = Pars::alpha_p*_u_star[i][j] + (1-Pars::alpha_p)*u[i][j];
				// err_u   += abs(_dt/_dx*(_p_n[i+1][j] - _p_n[i][j]));
				err_u   += abs(Pars::rho*_dy*_dx*(u[i+1][j] - u[i][j]));
			}
		}
		// printf("\n%f\n", err_u);

		for (int i = 0; i < _nx; i++) // excluding boundary
		{
			for (int j = 0; j < _ny-2; j++) // excluding boundary
			{
				// v[i][j] += -_dt/_dy*(_p_n[i][j] - _p_n[i][j-1]);
				// err_v   += abs(_dt/_dy*(_p_n[i][j] - _p_n[i][j-1]));
				// v[i][j] = _v_star[i][j] - _dt/(_dy*Pars::rho)*(_p_n(i,j+1) - _p_n(i,j));
				// v[i][j] = Pars::alpha_p*_v_star[i][j] + (1-Pars::alpha_p)*v[i][j];
				// err_v   += abs(_dt/_dy*(_p_n[i][j] - _p_n[i][j-1]));
				err_v   += abs(Pars::rho*_dy*_dx*(v[i][j+1] - v[i][j]));
			}
		}

		// printf("\n%f %f %f \n", err_p, err_u, err_v);
		cout <<  "\radded mass: " << err_p << ", no. iterations: " << n_iter << ", error u: " << err_u << ", error v: " << err_v/*<< "\r"*/;
	} while ( err_p > 1.0e-3 /*|| err_u > 1.0e-4 || err_v > 1.0e-4*/ && n_iter < 100);

		u = _u_star;
		v = _v_star;

	printf("\ncalculation done {(-,-)}.\n");


	// == saving data ==
	std::ofstream outs;
	outs.open("output/data_p.dat");
	for (int i = 0; i < _nx; i++)
	{
		for (int j = 0; j < _ny; j++)
		{
			outs << x_IJ[i][j] << " " << y_IJ[i][j] << " " << p[i][j] /*_p_n(i,j)*/ << "\n";
		}
	}
	outs.close();

	outs.open("output/data_u.dat");
	for (int i = 0; i < _nx-1; i++)
	{
		for (int j = 0; j < _ny; j++)
		{
			outs << x_ij_u[i][j] << " " << y_ij_u[i][j] << " " << u[i][j] << "\n";
		}
	}
	outs.close();

	outs.open("output/data_v.dat");
	for (int i = 0; i < _nx; i++)
	{
		for (int j = 0; j < _ny-1; j++)
		{
			// outs << x_ij_u[i][j] << " " << y_ij_u[i][j] << " " << u[i][j] << "\n";
			outs << x_ij_v[i][j] << " " << y_ij_v[i][j] << " " << v[i][j] << "\n";
		}
	}
	outs.close();

}