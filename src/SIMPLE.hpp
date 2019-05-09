#ifndef INCLUDED_SIMPLE
#define INCLUDED_SIMPLE

#ifndef INCLUDED_GLOBAL
#include "../global.hpp"
#endif

// #include <locale> // std::toupper, std::tolower

class SIMPLE
{

public:

  void MomentumEq(const std::vector<std::vector<double>> p,const std::vector<std::vector<double>> u,
    const std::vector<std::vector<double>> v,std::vector<double> ul,std::vector<double> ur,
    std::vector<double> vt,std::vector<double> vb,const std::vector<int> idx_in,
    const std::vector<int> idx_out,std::vector<std::vector<double>> &u_star,
    std::vector<std::vector<double>> &v_star);



};


#endif
