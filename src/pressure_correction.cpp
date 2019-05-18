#include "SIMPLE.hpp"

// =====================================================
// =====================================================
void SIMPLE::pressure_correction(std::vector<std::vector<double>> &p)
{
	// == reset class variables
	// this->_p_n.resize(Pars::nx, std::vector<double>(Pars::ny, 0.0e0));
	// this->_p_n_1.resize(Pars::nx, std::vector<double>(Pars::ny, 0.0e0));

	// _p_n   = arma::ones(Pars::nx, Pars::ny)*1.0e-4;
	// _p_n_1 = arma::ones(Pars::nx, Pars::ny)*1.0e-4;

	// _p_n.col(0)          = zeros<vec>(Pars::ny);
	// _p_n.col(Pars::nx-1) = zeros<vec>(Pars::ny);
	// _p_n.row(0)          = zeros<rowvec>(Pars::nx);
	// _p_n.row(Pars::ny-1) = zeros<rowvec>(Pars::nx);

	// _p_n_1.col(0)          = zeros<vec>(Pars::ny);
	// _p_n_1.col(Pars::nx-1) = zeros<vec>(Pars::ny);
	// _p_n_1.row(0)          = zeros<rowvec>(Pars::nx);
	// _p_n_1.row(Pars::ny-1) = zeros<rowvec>(Pars::nx);

	// == internal variables
	double a,b,c,d;
	// double p_n_bar;
	double dphi_x(0), dphi_y(0);
	int n_it_x=0;
	int n_it_y=0;
	// -- independent value of internal variables
	a = 2*(_dt/pow(_dx,2) + _dt/pow(_dy,2));
	b = -_dt/pow(_dx,2);
	c = -_dt/pow(_dy,2);

	// // == calculating corrected pressure 
	// do
	// {
	// 	for (int i = 1; i < _nx-1; i++) // excluding boundary
	// 	{
	// 		for (int j = 1; j < _ny-1; j++) // excluding boundary
	// 		{
	// 			// -- dependent value of internal variables
	// 			d = Pars::rho/_dx*(_u_star[i][j] - _u_star[i-1][j]) + Pars::rho/_dy*(_v_star[i][j] - _v_star[i][j-1]);
	// 			p_n_bar = -1/a*(b*_p_n[i+1][j] + b*_p_n[i-1][j] + c*_p_n[i][j+1] + c*_p_n[i][j-1] + d);
	// 			_p_n_1[i][j] = _p_n[i][j] + Pars::omega*(p_n_bar - _p_n[i][j]);
	// 			// printf("%f\n", _p_n_1[i][j]);
	// 		}
	// 	}
	// 	n_iter += 1;
	// 	get_error(dphi);
 //    cout <<  "\r" << dphi << " " << n_iter/*<< "\r"*/;
 //    // fflush ( stdin );
	// 	// printf("%f\n", dphi);
	// 	// -- update pressure correction
	// 	_p_n = _p_n_1;
	// } while(dphi > 1.0e-4);


	// == horizontal swipe ==
	arma::mat tdma_x(_nx, _nx, arma::fill::zeros);
	arma::vec rhs_x(_nx, arma::fill::zeros);
	arma::vec sol_x(_nx);

	// -- construct lhs matrix
	for (int i = 0; i < _nx; i++)
	{
		if (i==0 || i==_nx-1)
		{
			tdma_x(i,i) = 1;
		}
		else
		{
			tdma_x(i,i)   = a;
			tdma_x(i,i-1) = b;
			tdma_x(i,i+1) = b;			
		}
	}
	// tdma_x.print();

	do
	{
		for (int j = 1; j < _ny-1; j++)
		{
			for (int i = 1; i < _nx-1; i++)
			{
				if (i==0 || i==_nx-1)
				{
					rhs_x(i)    = 0.0e0;
				}
				else
				{
					d = (Pars::rho/_dx)*(_u_star[i][j] - _u_star[i-1][j]) + (Pars::rho/_dy)*(_v_star[i][j] - _v_star[i][j-1]);
					rhs_x(i) = -1*(c*_p_n_1(i,j-1) + c*_p_n(i,j+1) + d);
				}

			}
			sol_x = solve(tdma_x, rhs_x);

			_p_n_1.col(j) = sol_x/*.t()*/; // !!! this is row !!!
		}
		// _p_n_1.row(1).print();

		n_it_x+=1;
		get_error(dphi_x);
    // cout <<  "\rerror swap x: " << dphi_x << ", no.iterations: " << n_it_x;
		// -- update value
		_p_n = _p_n_1;
	} while(dphi_x > 1.0e-4 && n_it_x < 500);

	// printf("\n");

	// == vertical swipe ==
	arma::mat tdma_y(_ny, _ny, arma::fill::zeros);
	arma::vec rhs_y(_ny, arma::fill::zeros);
	arma::vec sol_y(_ny);

	// -- construct lhs matrix
	for (int i = 0; i < _ny; i++)
	{
		if (i==0 || i==_ny-1)
		{
			tdma_y(i,i) = 1;
		}
		else
		{
			tdma_y(i,i)   = a;
			tdma_y(i,i-1) = c;
			tdma_y(i,i+1) = c;			
		}
	}
	// tdma_y.print();

	do
	{
		for (int i = 1; i < _nx-1; i++)
		{
			for (int j = 0; j < _ny; j++)
			{
				if (j==0 || j==_ny-1)
				{
					rhs_y(j)    = 0.0e0;
				}
				else
				{
					d = Pars::rho/_dx*(_u_star[i][j] - _u_star[i-1][j]) + Pars::rho/_dy*(_v_star[i][j] - _v_star[i][j-1]);
					rhs_y(j) = -(b*_p_n_1(i-1,j) + b*_p_n(i+1,j) + d);
				}
			}
			// tdma_x.print();
			// cout << tdma_x.n_cols;
			sol_y = solve(tdma_y, rhs_y);

			// sol_y.print();
			// printf("\n");
			_p_n_1.row(i) = sol_y.t();
		}

		n_it_y+=1;
		get_error(dphi_y);
    // cout <<  "\rerror swap y: " << dphi_y << ", no.iterations: " << n_it_y;
		// -- update value
		_p_n = _p_n_1;
	} while(dphi_y > 1.0e-4 && n_it_y < 500);


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
			 // dphi += std::abs(_p_n_1[i][j] - _p_n[i][j]);
			 // sum_dp += std::abs(_p_n_1[i][j]);
 			 dphi   += std::abs(_p_n_1(i,j) - _p_n(i,j));
			 sum_dp += std::abs(_p_n_1(i,j));
		}
	}	
	dphi /= sum_dp;

}
// =====================================================
// =====================================================