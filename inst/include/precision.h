
#ifndef ___PRECISION___
#define ___PRECISION___

#include <RcppArmadillo.h>

namespace ostats
{

struct precision
{
  arma::sp_mat Q;
  std::vector<std::vector<int>> nbs;
  std::vector<std::vector<int>> extended_nbs;
  precision(const arma::sp_mat&,
            const std::vector<std::vector<int>>&,
            const std::vector<std::vector<int>>&);
};


std::vector<std::vector<int>> lower_nbs(const arma::sp_mat&);
std::vector<std::vector<int>> extended_lower_nbs(std::vector<std::vector<int>>);

}

#endif // ___PRECISION___
