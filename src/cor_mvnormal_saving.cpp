
# include "cor_mvnormal_saving.h"

using namespace bqp;
using namespace Rcpp;


// namespace
namespace ostats
{

  double dense_mvnormal_savings(const arma::mat& mean_x, const arma::sp_mat& Q, const double& n) {
    return n * arma::as_scalar(mean_x.t() * Q * mean_x);
  }

  arma::sp_mat BQP_A(const double& n, const arma::mat& mean_x, const arma::sp_mat& Q) {
    return - n * (mean_x * mean_x.t()) % Q;
  }

  arma::mat BQP_b(const double& n, const arma::mat& mean_x, const arma::sp_mat& Q, const double& beta) {
    return 2.0 * n * mean_x % (Q * mean_x) - beta;
  }

  Rcpp::List format_results(const BQP_res& res) {
      return Rcpp::List::create(_["S_max"] = res.max_value,
                                _["J_max"] = res.max_subset,
                                _["k_min"] = res.min_subset_size);
  }

} // unnamed namespace

namespace ostats
{

  cor_mvnormal_saving_state::cor_mvnormal_saving_state
  (
    const int& w,
    const int& p,
    precision&& precision_,
    linear_const_penalty& penalty_
  ) :
    vec_moments_state(w, p),
    prec_obj(precision_),
    penalty(penalty_) {}


  bqp::BQP_res cor_mvnormal_saving(const arma::mat& mean,
                                   const int& m,
                                   const precision& precision,
                                   linear_const_penalty& penalty)
	{
    int p = mean.n_elem;
    bqp::BQP_res sparse_res;
    bqp::BQP_res dense_res;

    if (penalty.beta > 0.0)
    {
      bqp::BQP curr_BQP(BQP_A(m, mean, precision.Q),
                        BQP_b(m, mean, precision.Q, penalty.beta),
                        - penalty.alpha_linear,
                        penalty.k_star(),
                        precision.nbs,
                        precision.extended_nbs);
      sparse_res = bqp::BQP_optim(curr_BQP);
    }

    double inf = std::numeric_limits<double>::infinity();
    if (penalty.alpha_const < inf)
    {
      double S_dense = dense_mvnormal_savings(mean, precision.Q, m) - penalty.alpha_const;
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

}  // namespace ostats

// [[Rcpp::export]]
Rcpp::List optimise_mvnormal_saving(const arma::mat& x,
                                    const arma::sp_mat& Q,
                                    const std::vector<std::vector<int>>& nbs,
                                    const std::vector<std::vector<int>>& extended_nbs,
                                    const double& alpha_dense = 0.0,
                                    const double& beta = 0.0,
                                    const double& alpha_sparse = 0.0) {
  arma::mat mean_x = mean(x).t();  // Column vector.
  double n = x.n_rows;
  ostats::precision prec_obj(Q, nbs, extended_nbs);
  ostats::linear_const_penalty penalty(alpha_dense, alpha_sparse, beta);
  return(ostats::format_results(cor_mvnormal_saving(mean_x, n, prec_obj, penalty)));
}  // optimise_mvnormal_savings

// [[Rcpp::export]]
double dense_mvnormal_savings(const arma::mat& mean_x, const arma::sp_mat& Q, const double& n) {
  return n * arma::as_scalar(mean_x.t() * Q * mean_x);
}






// // [[Rcpp::export]]
// Rcpp::List optimise_mvnormal_saving_old(const arma::mat& x,
//                                         const arma::sp_mat& Q,
//                                         const std::vector<std::vector<int>>& nbs,
//                                         const std::vector<std::vector<int>>& extended_nbs,
//                                         const double& alpha_dense = 0.0,
//                                         const double& beta = 0.0,
//                                         const double& alpha_sparse = 0.0) {
//   Rcpp::List sparse_res;
//   Rcpp::List dense_res;
//   double k_star;
//   arma::mat mean_x = mean(x).t();  // Column vector.
//
//   if (beta > 0.0)
//   {
//     double n = x.n_rows;
//     k_star = BQP_k_star(alpha_dense, alpha_sparse, beta);
//     BQP curr_BQP(BQP_A(n, mean_x, Q),
//                  BQP_b(n, mean_x, Q, beta),
//                  - alpha_sparse,
//                  k_star,
//                  nbs,
//                  extended_nbs);
//     sparse_res = format_results(BQP_optim(curr_BQP));
//   }
//   if (alpha_dense < R_PosInf)
//   {
//     double S_dense = dense_mvnormal_savings(mean_x, Q, x.n_rows) - alpha_dense;
//     dense_res = Rcpp::List::create(_["S_max"] = S_dense,
//                                    _["J_max"] = seq_int(1, x.n_cols));
//   }
//
//   if (beta > 0.0 && Rcpp::traits::is_infinite<REALSXP>(alpha_dense))
//   {
//     return sparse_res;
//   }
//   else if (beta <= 0.0 && alpha_dense < R_PosInf)
//   {
//     return dense_res;
//   }
//   else if (beta > 0.0 && alpha_dense < R_PosInf)
//   {
//     double k_min = sparse_res["k_min"];
//     if (k_min > k_star) return dense_res;
//     else
//     {
//       double sparse_S_max = sparse_res["S_max"];
//       double dense_S_max = dense_res["S_max"];
//       if (sparse_S_max > dense_S_max) return sparse_res;
//       else return dense_res;
//     }
//   }
//   else
//   {
//     Rcpp::List empty_list = List::create();
//     return empty_list;
//   }
// }  // optimise_mvnormal_savings_old



