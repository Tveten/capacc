
#ifndef ___PENALTY___
#define ___PENALTY___

#include <cmath>
#include <iostream>
#include <limits>

namespace ostats
{

  struct penalty
  {
    virtual double P(const int& k) = 0;
  };

  struct constant_penalty : virtual public penalty
  {
    double alpha_const;
    constant_penalty(const double&);
    constant_penalty(const constant_penalty&);
    double P(const int&);
  };

  struct linear_penalty : virtual public penalty
  {
    double alpha_linear;
    double beta;
    linear_penalty(const double&, const double&);
    linear_penalty(const linear_penalty&);
    double P(const int&);
  };

  struct linear_const_penalty : public constant_penalty, public linear_penalty
  {
    linear_const_penalty(const double&, const double&, const double&);
    linear_const_penalty(const linear_const_penalty&);
    double k_star();
    double P(const int&);
  };

  linear_const_penalty mvnormal_default_penalty(const int&, const int&, const double& b = 1.0);
  linear_const_penalty mvnormal_default_point_penalty(const int&, const int&, const double& b = 1.0);

}

#endif // ___PENALTY___
