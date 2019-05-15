#include "SIMPLE.hpp"

// =====================================================
// =====================================================

SIMPLE::SIMPLE()
{
	this->_nx = Pars::nx;
	this->_ny = Pars::ny;
	this->_dx = Pars::dx;
	this->_dy = Pars::dy;
	this->_dt = Pars::dt;
	// == initialization of internal variables for pressure correction purpose
	// this->_p_n.resize(Pars::nx, std::vector<double>(Pars::ny, 0.0e0));
	// this->_p_n_1.resize(Pars::nx, std::vector<double>(Pars::ny, 0.0e0));
}

// =====================================================
// =====================================================

void SIMPLE::MomentumEq(const std::vector<std::vector<double>> p,const std::vector<std::vector<double>> u,
    const std::vector<std::vector<double>> v, std::vector<double> ul, std::vector<double> ur,
    std::vector<double> vt, std::vector<double> vb,const std::vector<int> idx_in,
    const std::vector<int> idx_out)
 {
    double dt = Pars::dt;

    // Ghost Points
	  ul.resize(Pars::ny);ur.resize(Pars::ny); // to assist x-momentum solving 
	  vt.resize(Pars::nx);vb.resize(Pars::nx); // to assist y-momentum solving 
 
    // Left and Right Boundary
    for(size_t j=0;j<Pars::ny;j++){
      if(idx_in[j]==1){
        ul[j]=Pars::u_i;
      }else{
        ul[j]=u[0][j];
      }
      if(idx_out[j]==1){ 
        // ur[j]=u[Pars::nx-2][j]; //? Velocity at Outlet <-- less safe
        ur[j] = 2.0*u[Pars::nx-2][j]-u[Pars::nx-3][j]; //Linear Extrapolation
      }else{
        ur[j]=u[Pars::nx-2][j];
      } 
    }
    // Top and Bottom Boundary 
    for(size_t i;i<Pars::nx;i++){
      vt[i] = v[i][Pars::nx-1];
      vb[i] = v[i][0];
    }

	// Momentum Equation x
	this->_u_star = u;  // u^{n+1}
	for(size_t i=0;i<Pars::nx-1;i++){
		for(size_t j=1;j<Pars::ny-1;j++){
		
		double vbar_1,vbar_2;
		vbar_1 = 0.5*(v[i][j]+v[i+1][j]);		// vbar_{i}{j+1/2}
		vbar_2 = 0.5*(v[i][j-1]+v[i+1][j-1]);	// vbar_{i}{j-1/2}
		
		double A,A1,A2,A3,A4;
		if(i==0){	// Points Near Left Boundary
			A1 = (pow(u[i+1][j],2)-pow(ul[j],2))/(2*Pars::dx);
			A2 = (u[i][j+1]*vbar_1-u[i][j-1]*vbar_2)/(2*Pars::dy);
			A3 = (u[i+1][j]-2*u[i][j]+ul[j])/pow(Pars::dx,2);
			A4 = (u[i][j+1]-2*u[i][j]+u[i][j-1])/pow(Pars::dy,2);			
		}
		else if(i==Pars::nx-2){ // Points Near Right Boundary
			A1 = (pow(ur[j],2)-pow(u[i-1][j],2))/(2*Pars::dx);
			A2 = (u[i][j+1]*vbar_1-u[i][j-1]*vbar_2)/(2*Pars::dy);
			A3 = (ur[j]-2*u[i][j]+u[i-1][j])/pow(Pars::dx,2);
			A4 = (u[i][j+1]-2*u[i][j]+u[i][j-1])/pow(Pars::dy,2);
		}else{
			A1 = (pow(u[i+1][j],2)-pow(u[i-1][j],2))/(2*Pars::dx);
			A2 = (u[i][j+1]*vbar_1-u[i][j-1]*vbar_2)/(2*Pars::dy);
			A3 = (u[i+1][j]-2*u[i][j]+u[i-1][j])/pow(Pars::dx,2);
			A4 = (u[i][j+1]-2*u[i][j]+u[i][j-1])/pow(Pars::dy,2);
		}
	
		A = -Pars::rho*(A1+A2)+Pars::miu*(A3+A4);

		_u_star[i][j] = Pars::rho*u[i][j]+ A*dt-((dt/Pars::dx)*(p[i+1][j]-p[i][j]));
		_u_star[i][j] = (1.0/Pars::rho)*_u_star[i][j];

		}
	}	

	// Momentum Equation y
	this->_v_star = v;  // v^{n+1}
	for(size_t i=1;i<Pars::nx-1;i++){
		for(size_t j=0;j<Pars::ny-1;j++){
		 
		 double ubar_1,ubar_2;
		 double B,B1,B2,B3,B4;

		//  if(i!=Pars::nx-1){
			
			ubar_1 = 0.5*(u[i-1][j]+u[i-1][j+1]);	// ubar_{i-1/2}{j}
			ubar_2 = 0.5*(u[i][j]+u[i][j+1]);		// ubar_{i+1/2}{j}
			
			if(j==0){	// Points Near Bottom Boundary
				B1 = (v[i+1][j]*ubar_2-v[i-1][j]*ubar_1)/(2*Pars::dx);
				B2 = (pow(v[i][j+1],2)-pow(vb[i],2))/(2*Pars::dy);
				B3 = (v[i+1][j]-2*v[i][j]+v[i-1][j])/pow(Pars::dx,2);
				B4 = (v[i][j+1]-2*v[i][j]+vb[i])/pow(Pars::dy,2);
			}else if(j==Pars::ny-2){	// Points Near Top Boundary
				B1 = (v[i+1][j]*ubar_2-v[i-1][j]*ubar_1)/(2*Pars::dx);
				B2 = (pow(vt[i],2)-pow(v[i][j-1],2))/(2*Pars::dy);
				B3 = (v[i+1][j]-2*v[i][j]+v[i-1][j])/pow(Pars::dx,2);
				B4 = (vt[i]-2*v[i][j]+v[i][j-1])/pow(Pars::dy,2);
			}else{
				B1 = (v[i+1][j]*ubar_2-v[i-1][j]*ubar_1)/(2*Pars::dx);
				B2 = (pow(v[i][j+1],2)-pow(v[i][j-1],2))/(2*Pars::dy);
				B3 = (v[i+1][j]-2*v[i][j]+v[i-1][j])/pow(Pars::dx,2);
				B4 = (v[i][j+1]-2*v[i][j]+v[i][j-1])/pow(Pars::dy,2);
			}
		//  }
	 //  else{	// Left Boundary (Outlet)
			
		// 	if(idx_out[j]==1){ //Outlet        
    //     ubar_1 = 0.5*(u[i-1][j]+u[i-1][j+1]);	// ubar_{i-1/2}{j}
		// 		ubar_2 = 0.5*(ur[j]+ur[j+1]);			// ubar_{i+1/2}{j}
		// 		if(j==0){	
		// 			B1 = (v[i][j]*ubar_2-v[i-1][j]*ubar_1)/(2*Pars::dx); //? Assume : v[nx][j]==v[nx-1][j]
		// 			B2 = (pow(v[i][j+1],2)-pow(vb[i],2))/(2*Pars::dy);
		// 			B3 = (v[i][j]-2*v[i][j]+v[i-1][j])/pow(Pars::dx,2);  //? Assume : v[nx][j]==v[nx-1][j]
		// 			B4 = (v[i][j+1]-2*v[i][j]+vb[i])/pow(Pars::dy,2);
		// 		}else{	
		// 			B1 = (v[i][j]*ubar_2-v[i-1][j]*ubar_1)/(2*Pars::dx); //? Assume : v[nx][j]==v[nx-1][j]
		// 			B2 = (pow(v[i][j+1],2)-pow(v[i][j-1],2))/(2*Pars::dy);
		// 			B3 = (v[i][j]-2*v[i][j]+v[i-1][j])/pow(Pars::dx,2);  //? Assume : v[nx][j]==v[nx-1][j]
		// 			B4 = (v[i][j+1]-2*v[i][j]+v[i][j-1])/pow(Pars::dy,2);
		// 		}
		// 	}
	 //  }

			B = -Pars::rho*(B1+B2)+Pars::miu*(B3+B4);

			_v_star[i][j] = Pars::rho*v[i][j]+B*dt-((dt/Pars::dy)*(p[i][j+1]-p[i][j]));
			_v_star[i][j] = (1.0/Pars::rho)*_v_star[i][j];
		}
	}
  // Left Boundary (Outlet) <-- does it means right boundary? since it is indexed from (nx-1)
  //Linear Extrapolation
	  for(size_t j=0;j<Pars::ny-1;j++){
	    if(idx_out[j]==1){
	      _v_star[Pars::nx-1][j] = 2.0*_v_star[Pars::nx-2][j]-_v_star[Pars::nx-3][j];
	    }
	  }

 }
