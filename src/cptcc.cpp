
#include "cptcc.h"

using namespace bqp;
using namespace Rcpp;

// namespace
namespace
{

  double dense_mvnormal_lr(const double& n1,
                           const double& n2,
                           const arma::mat& mean1,
                           const arma::mat& mean2,
                           const arma::sp_mat& Q) {
    return n1 * arma::as_scalar(mean1.t() * Q * mean1) +
             n2 * arma::as_scalar(mean2.t() * Q * mean2);
  }

  arma::sp_mat BQP_A(const double& n1,
                     const double& n2,
                     const arma::mat& mean1,
                     const arma::mat& mean2,
                     const arma::sp_mat& Q) {
    return - n1 * (mean1 * mean1.t()) % Q - n2 * (mean2 * mean2.t()) % Q;
  }

  arma::mat BQP_B(const double& n1,
                  const double& n2,
                  const arma::mat& mean1,
                  const arma::mat& mean2,
                  const arma::sp_mat& Q,
                  const double& beta) {
    return 2.0 * (n1 * mean1 % (Q * mean1) + n2 * mean2 % (Q * mean2))  - beta;
  }

  std::vector<int> seq_int(const int& start, const int& end)
  {
    std::vector<int> v(end - start + 1);
    std::iota(v.begin(), v.end(), start);
    return v;
  }

  arma::mat centralise(arma::mat x)
  {
    arma::mat mean = arma::mean(x, 0);
    x.each_row() -= mean;
    return x;
  }

  Rcpp::List format_cpt_results(const BQP_res& res)
  {
      return Rcpp::List::create(_["S_max"] = res.max_value,
                                _["J_max"] = res.max_subset);
  }
} // unnamed namespace

namespace ostats {
  bqp::BQP_res cor_mvnormal_lr(const arma::mat& mean1,
                               const arma::mat& mean2,
                               const int& m1,
                               const int& m2,
                               const precision& precision,
                               linear_const_penalty& penalty)
  {
    int p = mean1.n_elem;
    bqp::BQP_res sparse_res;
    bqp::BQP_res dense_res;

    if (penalty.beta > 0.0)
    {
      bqp::BQP curr_BQP(BQP_A(m1, m2, mean1, mean2, precision.Q),
                        BQP_B(m1, m2, mean1, mean2, precision.Q, penalty.beta),
                        - penalty.alpha_linear,
                        penalty.k_star(),
                        precision.nbs,
                        precision.extended_nbs);
      sparse_res = bqp::BQP_optim(curr_BQP);
    }

    double inf = std::numeric_limits<double>::infinity();
    if (penalty.alpha_const <  inf)
    {
      double S_dense = dense_mvnormal_lr(m1, m2, mean1, mean2, precision.Q) - penalty.alpha_const;
      std::vector<int> J_dense(p);
      std::iota(J_dense.begin(), J_dense.end(), 1);
      dense_res = bqp::BQP_res(S_dense, J_dense, p);
    }

    if (penalty.beta > 0.0 && std::isinf(penalty.alpha_const)) return sparse_res;
    else if (penalty.beta <= 0.0 && penalty.alpha_const < inf) return dense_res;
    else if (penalty.beta > 0.0 && penalty.alpha_const < inf)
    {
      if (sparse_res.min_subset_size > penalty.k_star()) return dense_res;
      else
      {
        if (sparse_res.max_value > dense_res.max_value) return sparse_res;
        else return dense_res;
      }
    }
    else
    {
      bqp::BQP_res empty_BQP_res;
      return empty_BQP_res;
    }
  }

}

using namespace ostats;

// [[Rcpp::export]]
Rcpp::List optimise_mvnormal_lr(const int& cpt,
                                const arma::mat& x,
                                const arma::sp_mat& Q,
                                const double& b = 1)
{
  int n = x.n_rows;
  int p = x.n_cols;
  std::vector<std::vector<int>> nbs = lower_nbs(Q);
  precision precision_obj(Q, nbs, extended_lower_nbs(nbs));
  linear_const_penalty penalty = mvnormal_default_penalty(n, p, b);
  arma::mat x_centered = centralise(x);
  arma::mat M = cumsum(x_centered, 0).t();
  arma::mat mean1 = M.col(cpt - 1) / cpt;
  arma::mat mean2 = (M.col(n - 1) - M.col(cpt - 1)) / (n - cpt);
  bqp::BQP_res lr = cor_mvnormal_lr(mean1, mean2, cpt, n - cpt, precision_obj, penalty);
  return(format_cpt_results(lr));
}

// [[Rcpp::export]]
Rcpp::List cptcc(const arma::mat& x,
                 const arma::sp_mat& Q,
                 const double& b = 1,
                 const int& min_seg_len = 2)
{
  int n = x.n_rows;
  int p = x.n_cols;
  std::vector<std::vector<int>> nbs = lower_nbs(Q);
  precision precision_obj(Q, nbs, extended_lower_nbs(nbs));
  linear_const_penalty penalty = mvnormal_default_penalty(n, p, b);
  arma::mat x_centered = centralise(x);
  arma::mat M = cumsum(x_centered, 0).t();

  std::vector<double> C(n + 1);
  std::vector<std::vector<int>> J(n + 1);
  std::vector<int> T = seq_int(min_seg_len, n - min_seg_len);
  #ifdef _OPENMP
    std::vector<bqp::BQP_res> lr(n);
    int n_threads = std::thread::hardware_concurrency();
    if (n_threads > 16) n_threads = n_threads / 2;
    else if (n_threads > 8 && n_threads <= 16) n_threads = 8;
    #pragma omp parallel for num_threads(n_threads)
    for (int i = 0; i < T.size(); i++)
    {
      int t = T[i];
      arma::mat mean1 = M.col(t - 1) / t;
      arma::mat mean2 = (M.col(n - 1) - M.col(t - 1)) / (n - t);
      lr[t] = cor_mvnormal_lr(mean1, mean2, t, n - t, precision_obj, penalty);
  	  C[t] = lr[t].max_value;
  	  J[t] = lr[t].max_subset;
    }
  #else
    bqp::BQP_res lr;
    for (auto& t : T)
    {
      arma::mat mean1 = M.col(t - 1) / t;
      arma::mat mean2 = (M.col(n - 1) - M.col(t - 1)) / (n - t);
      lr = cor_mvnormal_lr(mean1, mean2, t, n - t, precision_obj, penalty);
  	  C[t] = lr.max_value;
  	  J[t] = lr.max_subset;
    }
  #endif

  double C_max = - std::numeric_limits<double>::infinity();
  int t_max;
  for (int t = T.front(); t <= T.back(); t++)
  {
    if (C[t] > C_max)
    {
      C_max = C[t];
      t_max = t;
    }
  }
  return Rcpp::List::create(Rcpp::_["cpt"] = t_max,
                            Rcpp::_["J"]   = J[t_max],
                            Rcpp::_["value"]  = C_max);
}

