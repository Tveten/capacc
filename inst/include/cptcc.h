
#include "bqp.h"
#include "precision.h"
#include "penalty.h"

#include <RcppArmadillo.h>
#include <limits>
#include <thread>
#ifdef _OPENMP
#include <omp.h>
#endif


namespace ostats
{
  bqp::BQP_res cor_mvnormal_lr(const arma::mat&,
                               const arma::mat&,
                               const int&,
                               const int&,
                               const precision&,
                               linear_const_penalty&);
}
