
#ifndef ___MVNORMAL_SAVING_H___
#define ___MVNORMAL_SAVING_H___

#include "moments.h"
#include "bqp.h"
#include "precision.h"
#include "penalty.h"

#include <RcppArmadillo.h>
#include <memory>
#include <list>
#include <limits>
// #include <iostream>
#include <thread>
#include <algorithm>
#include <boost/functional/hash.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace ostats
{

  struct cor_mvnormal_saving_state: vec_moments_state
  {
    precision prec_obj;
    linear_const_penalty penalty;
    cor_mvnormal_saving_state(const int&,
                              const int&,
                              // const precision&,
                              precision&&,
                              linear_const_penalty&);
  };

  bqp::BQP_res cor_mvnormal_saving(const arma::mat&,
                                   const int&,
                                   const precision&,
                                   linear_const_penalty&);

  template <class state_type>
    state_type cor_mvnormal_savings(state_type S)
    {
      auto n = S.M1.size() - 1;
      #ifdef _OPENMP
        std::vector<bqp::BQP_res> saving(n);
        int n_threads = std::thread::hardware_concurrency();
        // if (n_threads > 16) n_threads = n_threads / 2;
        // else if (n_threads > 8 && n_threads <= 16) n_threads = 8;
        if (n_threads > 8) n_threads = 8;
        #pragma omp parallel for num_threads(n_threads)
        for (int i = 0; i < S.T.size(); i++)
        {
          int t = S.T[i];
  	      arma::mat mean = (S.M1[n] - S.M1[t]) / (n - t);  // Mean of obs. t + 1, ..., n.
          saving[t] = cor_mvnormal_saving(mean, n - t, S.prec_obj, S.penalty);
      	  S.C[t] = - saving[t].max_value;  // Summing over obs (i + 1) to n.
      	  S.J[t] = saving[t].max_subset;
        }
      #else
        bqp::BQP_res saving;
        for (auto& t : S.T) {
  	      arma::mat mean = (S.M1[n] - S.M1[t]) / (n - t);  // Mean of obs. t + 1, ..., n.
          saving = cor_mvnormal_saving(mean, n - t, S.prec_obj, S.penalty);
      	  S.C[t] = - saving.max_value;  // Summing over obs (i + 1) to n.
      	  S.J[t] = saving.max_subset;
        }
      #endif
//
      return std::move(S);
    }  // cor_mvnormal_savings

  double dense_mvnormal_savings(const arma::mat& mean_x, const arma::sp_mat& Q, const double& n);

  arma::sp_mat BQP_A(const double& n, const arma::mat& mean_x, const arma::sp_mat& Q);

  arma::mat BQP_b(const double& n, const arma::mat& mean_x, const arma::sp_mat& Q, const double& beta);

} // namespace ostats

#endif // ___MVNORMAL_SAVING_H___

