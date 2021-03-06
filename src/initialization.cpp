#include "initialization.hpp"

// =====================================================
// =====================================================
Initialization::Initialization()
{
	this->_nx = Pars::nx;
  this->_ny = Pars::ny;
  this->_lx = Pars::lx;
  this->_ly = Pars::ly;
  this->_dx = Pars::dx;
  this->_dy = Pars::dy;
}
// =====================================================
// =====================================================
Initialization::~Initialization()
{

}
// =====================================================
// =====================================================
void Initialization::init_domain(std::vector<std::vector<double>> &x_IJ, std::vector<std::vector<double>> &y_IJ,
   std::vector<std::vector<double>> &x_ij_u, std::vector<std::vector<double>> &y_ij_u,
   std::vector<std::vector<double>> &x_ij_v, std::vector<std::vector<double>> &y_ij_v,
   std::vector<int> &idx_in,std::vector<int> &idx_out,
   std::vector<std::vector<double>> &P,std::vector<std::vector<double>>&u, std::vector<std::vector<double>>&v)
 {
   // -- grid coordinate system
   x_IJ.resize(_nx, std::vector<double>(_ny));
   y_IJ.resize(_nx, std::vector<double>(_ny));
   x_ij_u.resize(_nx-1, std::vector<double>(_ny));
   y_ij_u.resize(_nx-1, std::vector<double>(_ny));
   x_ij_v.resize(_nx, std::vector<double>(_ny-1));
   y_ij_v.resize(_nx, std::vector<double>(_ny-1));
   // -- boundary indexes
   idx_in.resize(_ny);
   idx_out.resize(_ny);
   // -- primitive variables
   P.resize(_nx, std::vector<double>(_ny));
   u.resize(_nx-1, std::vector<double>(_ny));
   v.resize(_nx, std::vector<double>(_ny-1));

   // == saving data ==
   std::ofstream outx;
   std::ofstream outy;
   std::ofstream outs;

   // -- setup for P domain
   outx.open("output/x_IJ.csv");
   outy.open("output/y_IJ.csv");
   outs.open("output/IJ.csv");
   for (size_t i = 0; i < _nx; i++) {
     for (size_t j = 0; j < _ny; j++) {
       x_IJ[i][j] = 0.0e0 + (double)i/((double)_nx-1)*_lx;
       y_IJ[i][j] = 0.0e0 + (double)j/((double)_ny-1)*_ly;
  		 outx << x_IJ[i][j] << ",";
       outy << y_IJ[i][j] << ",";
       outs << x_IJ[i][j] << "," << y_IJ[i][j] << "\n";
     }
     outx << "\n";
     outy << "\n";
   }
   outx.close();
   outy.close();
   outs.close();
  // -- setup for u domain, forward staggered
  outx.open("output/x_ij_u.csv");
  outy.open("output/y_ij_u.csv");
  outs.open("output/ij_u.csv");
   for (size_t i = 0; i < _nx-1; i++) {
     for (size_t j = 0; j < _ny; j++) {
       x_ij_u[i][j] = x_IJ[i][j] + 0.5*_dx;
       y_ij_u[i][j] = y_IJ[i][j];
       outx << x_ij_u[i][j] << ",";
       outy << y_ij_u[i][j] << ",";
       outs << x_ij_u[i][j] << "," << y_ij_u[i][j] << "\n";
     }
     outx << "\n";
     outy << "\n";
   }
   outx.close();
   outy.close();
   outs.close();
  // -- setup for v domain, forward staggered
  outx.open("output/x_ij_v.csv");
  outy.open("output/y_ij_v.csv");
  outs.open("output/ij_v.csv");
   for (size_t i = 0; i < _nx; i++) {
     for (size_t j = 0; j < _ny-1; j++) {
       x_ij_v[i][j] = x_IJ[i][j];
       y_ij_v[i][j] = y_IJ[i][j] + 0.5*_dy;
       outx << x_ij_v[i][j] << ",";
       outy << y_ij_v[i][j] << ",";
       outs << x_ij_v[i][j] << "," << y_ij_v[i][j] << "\n";
     }
     outx << "\n";
     outy << "\n";
   }
   outx.close();
   outy.close();
   outs.close();
 	// outs.open("output/IJ.csv");
 	// for (int i = 0; i < x.size(); i++)
 	// 	outs << x[i] << " " << y[i] << "\n";
 	// outs.close();

  // Inlet and Outlet 
  for(size_t i=0;i<_ny;i++){
    // -- 0.5 < y < 0.75 && x=0
    if((y_IJ[0][i]>=(_lx-Pars::l_i-0.25))&&(y_IJ[0][i]<=(_lx-0.25))){
      idx_in[i] = 1;
    }else{
      idx_in[i] = 0;
    }
    // -- 0 < y < 0.25 && x=1
    if(y_IJ[_ny-1][i]<=Pars::l_o){
      idx_out[i] = 1;
    }else{
      idx_out[i] = 0;
    }
  }

  // Initial Pressure and Velocity Value
    // Pressure
    for(size_t i=0;i<_nx;i++){
      for(size_t j=0;j<_ny;j++){
        // -- inlet
        if((i==0)&&(idx_in[j]==1)){
          P[i][j] = Pars::p_i;
        // -- outlet  
        }else if((x_IJ[i][j]==_lx)&&(y_IJ[i][j]<=Pars::l_o)){
          P[i][j] = Pars::p_o;
        }
        // -- otherwise (inner points, including wall). however, p_wall will never be used explicitly. hence it can be stated safe to make these values equal to zeros.
        else{
          P[i][j] = 0.0;
        }
      }
    }
    // Velocity (u)
    for(size_t i=0;i<_nx-1;i++){
      for(size_t j=0;j<_ny;j++){
        u[i][j] = 0.0;
      }
    }
    // Velocity (v)
    for(size_t i=0;i<_nx;i++){
      for(size_t j=0;j<_ny-1;j++){
        if((i==0)&&(y_ij_v[0][j]>=(_lx-Pars::l_i-0.25))&&(y_ij_v[0][j]<=(_lx-0.25))){
          v[i][j] = Pars::v_i;
        }
        else{
          v[i][j] = 0.0;
        }
      }
    }

 }
