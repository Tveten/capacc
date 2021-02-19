
#include "penalty.h"

namespace ostats
{
  constant_penalty::constant_penalty(const double& alpha_const_)
    : alpha_const(alpha_const_) {}

  constant_penalty::constant_penalty(const constant_penalty& penalty)
    : alpha_const(penalty.alpha_const) {}

  double constant_penalty::P(const int& k)
  {
    return alpha_const;
  }

  linear_penalty::linear_penalty(const double& alpha_linear_, const double& beta_)
    : alpha_linear(alpha_linear_), beta(beta_) {}

  linear_penalty::linear_penalty(const linear_penalty& penalty)
    : alpha_linear(penalty.alpha_linear), beta(penalty.beta) {}

  double linear_penalty::P(const int& k)
  {
    return alpha_linear + beta * k;
  }

  linear_const_penalty::linear_const_penalty
  (
    const double& alpha_const_,
    const double& alpha_linear_,
    const double& beta_
  ) :
    constant_penalty(alpha_const_),
    linear_penalty(alpha_linear_, beta_) {}

  linear_const_penalty::linear_const_penalty
  (
    const linear_const_penalty& penalty
  ) :
    constant_penalty(penalty.alpha_const),
    linear_penalty(penalty.alpha_linear, penalty.beta) {}

  double linear_const_penalty::k_star()
  {
    if (beta > 0.0) return (alpha_const - alpha_linear) / beta;
    else return 0.0;
  }

  double linear_const_penalty::P(const int& k)
  {
    if (k < k_star()) return alpha_linear + beta * k;
    else return alpha_const;
  }

  linear_const_penalty mvnormal_default_penalty(const int& n, const int& p, const double& b)
  {
    double psi = 2 * log(n);
    double a_const = b * (p + 2 * psi + 2 * sqrt(p * psi));
    double a_lin = b * 2 * psi;
    double beta = b * 2 * log(p);
    linear_const_penalty penalty(a_const, a_lin, beta);
    return penalty;
  }

  linear_const_penalty mvnormal_default_point_penalty(const int& n, const int& p, const double& b)
  {
    linear_const_penalty collective_penalty = mvnormal_default_penalty(n, p, b);
    double a_const = std::numeric_limits<double>::infinity();
    double a_lin = 0.0;
    double beta = collective_penalty.beta + collective_penalty.alpha_linear;
    linear_const_penalty penalty(a_const, a_lin, beta);
    return penalty;
  }

  linear_const_penalty mvnormal_default_sparse_penalty(const int& n, const int& p, const double& b)
  {
    double psi = 2 * log(n);
    double a_const = 100000000;
    double a_lin = b * 2 * psi;
    double beta = b * 2 * log(p);
    linear_const_penalty penalty(a_const, a_lin, beta);
    return penalty;
  }


} // namespace ostats

