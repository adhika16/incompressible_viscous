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
	
	double dt = 0.01;

	// Ghost Points
	ul.resize(Pars::ny);ur.resize(Pars::ny);
	vt.resize(Pars::nx);vb.resize(Pars::nx);
    // Left and Right Boundary
    for(size_t j=0;j<Pars::ny;j++){
      if(idx_in[j]==1){
        ul[j]=Pars::u_i;
      }else{
        ul[j]=0.0e0;
      }
      if(idx_out[j]==1){ 
        ur[j]=u[Pars::nx-2][j]; //? Velocity at Outlet
      }else{
        ur[j]=0.0e0;
      } 
    }

    // Top and Bottom Boundary 
    for(size_t i;i<Pars::nx;i++){
      vt[i] = 0.0e0;
      vb[i] = 0.0e0;
    }

	// Momentum Equation x
	std::vector<std::vector<double>> u_new = u;  // u^{n+1}
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

		u_new[i][j] = Pars::rho*u[i][j]+ A*dt-((dt/Pars::dx)*(p[i+1][j]-p[i][j]));
		u_new[i][j] = (1.0/Pars::rho)*u_new[i][j];

		}
	}	

	// Momentum Equation y
	std::vector<std::vector<double>> v_new = v;  // u^{n+1}
	for(size_t i=1;i<Pars::nx-1;i++){
		for(size_t j=0;j<Pars::ny-1;j++){
		 
		 double ubar_1,ubar_2;
		 double B,B1,B2,B3,B4;

		 if(i!=Pars::nx-1){
			
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
		 }
		 else{	// Left Boundary (Outlet)
			
			if(idx_out[j]==1){ //Outlet
				ubar_1 = 0.5*(u[i-1][j]+u[i-1][j+1]);	// ubar_{i-1/2}{j}
				ubar_2 = 0.5*(ur[j]+ur[j+1]);			// ubar_{i+1/2}{j}
				if(j==0){	
					B1 = (v[i][j]*ubar_2-v[i-1][j]*ubar_1)/(2*Pars::dx); //? Assume : v[nx][j]==v[nx-1][j]
					B2 = (pow(v[i][j+1],2)-pow(vb[i],2))/(2*Pars::dy);
					B3 = (v[i][j]-2*v[i][j]+v[i-1][j])/pow(Pars::dx,2);  //? Assume : v[nx][j]==v[nx-1][j]
					B4 = (v[i][j+1]-2*v[i][j]+vb[i])/pow(Pars::dy,2);
				}else{	
					B1 = (v[i][j]*ubar_2-v[i-1][j]*ubar_1)/(2*Pars::dx); //? Assume : v[nx][j]==v[nx-1][j]
					B2 = (pow(v[i][j+1],2)-pow(v[i][j-1],2))/(2*Pars::dy);
					B3 = (v[i][j]-2*v[i][j]+v[i-1][j])/pow(Pars::dx,2);  //? Assume : v[nx][j]==v[nx-1][j]
					B4 = (v[i][j+1]-2*v[i][j]+v[i][j-1])/pow(Pars::dy,2);
				}
			}
		 }

			B = -Pars::rho*(B1+B2)+Pars::miu*(B3+B4);

			v_new[i][j] = Pars::rho*v[i][j]+B*dt-((dt/Pars::dy)*(p[i][j+1]-p[i][j]));
			v_new[i][j] = (1.0/Pars::rho)*v_new[i][j];
		}
	}

	// == saving data ==
	// std::ofstream outs;
	// outs.open("output/airfoildata.dat");
	// for (int i = 0; i < x.size(); i++)
	// 	outs << x[i] << " " << y[i] << "\n";
	// outs.close();

	return EXIT_SUCCESS;
}
