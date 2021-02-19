#include "ostats.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <RcppArmadillo.h>

Rcpp::DataFrame format_capacc_output(std::list<std::tuple<int,int,int,double>> anoms)
{
  std::vector<int> start;
  std::vector<int> end;
  std::vector<int> variate;
  std::vector<double> size;
  for (auto& anom : anoms)
  {
    start.push_back(std::get<0>(anom));
    end.push_back(std::get<1>(anom));
    variate.push_back(std::get<2>(anom));
    size.push_back(std::get<3>(anom));
  }
  return Rcpp::DataFrame::create(Rcpp::_["start"] = start,
                                 Rcpp::_["end"] = end,
                                 Rcpp::_["variate"] = variate,
                                 Rcpp::_["size"] = size);
}

// [[Rcpp::export]]
Rcpp::DataFrame capacc(const arma::mat& x, const arma::sp_mat& Q,
                       const double& b = 1, const double& b_point = 1,
                       const int& min_seg_len = 2,
                       const int& max_seg_len = 100000000)
{

  // arma::mat x;
  //
  // // 100 x 3 - dimensional iid normal data with a large collective anomaly from time 50 to 69
  // // of size 10 in the first component only. Also a point anomaly at time 25 of size 50 in the
  // // first component.
  // x.load("./mv_testdata.csv");
  // std::cout << "Col1: " << x(0, 0) << " " << x(1, 0) << " " << x(2, 0) << " " <<std::endl;

  // initialise state
  double n = x.n_rows;
  int p = x.n_cols;
  std::vector<std::vector<int>> nbs = lower_nbs(Q);
  precision precision_obj(Q, nbs, extended_lower_nbs(nbs));
  linear_const_penalty penalty = mvnormal_default_penalty(n, p, b);
  // double psi = 2 * log(n);
  // double a_const = b * (p + 2 * psi + 2 * sqrt(p * psi));
  // linear_const_penalty penalty(a_const, 0.0, 0.0);
  linear_const_penalty point_penalty = mvnormal_default_point_penalty(n, p, b_point);
  cor_mvnormal_capa_state S(n, p, std::move(precision_obj), penalty, point_penalty,
                            min_seg_len, max_seg_len);

  arma::mat x_t = x.t();  // Functions below expect column vectors.

  // run capa on all the data
  for (int i = 0; i < n; i++)
  {
    S = capa(
          capa_cor_mvnormal(
            update_cpts(
              prune(
                op(
                  cor_mvnormal_savings(
                    update_candidate_cpts(
                      update_vec_moments(
                        update_observation(std::move(S), x_t.col(i))
    ))))))));
  }


  // show the results
  auto locs_and_vars = mvcapa_locations_and_variables(S);
  // std::cout << "***************" << std::endl;
  // for(auto& lv : locs_and_vars)
  //   {
  //     std::cout << "(" << std::get<0>(lv)  << "," << std::get<1>(lv) << "," << std::get<2>(lv) << ")" << std::endl;
  //   }

  return format_capacc_output(locs_and_vars);
}

// [[Rcpp::export]]
Rcpp::DataFrame capacc_sparse(const arma::mat& x, const arma::sp_mat& Q,
                              const double& b = 1, const double& b_point = 1,
                              const int& min_seg_len = 2,
                              const int& max_seg_len = 100000000)
{

  // arma::mat x;
  //
  // // 100 x 3 - dimensional iid normal data with a large collective anomaly from time 50 to 69
  // // of size 10 in the first component only. Also a point anomaly at time 25 of size 50 in the
  // // first component.
  // x.load("./mv_testdata.csv");
  // std::cout << "Col1: " << x(0, 0) << " " << x(1, 0) << " " << x(2, 0) << " " <<std::endl;

  // initialise state
  double n = x.n_rows;
  int p = x.n_cols;
  std::vector<std::vector<int>> nbs = lower_nbs(Q);
  precision precision_obj(Q, nbs, extended_lower_nbs(nbs));
  linear_const_penalty penalty = mvnormal_default_sparse_penalty(n, p, b);
  // double psi = 2 * log(n);
  // double a_const = b * (p + 2 * psi + 2 * sqrt(p * psi));
  // linear_const_penalty penalty(a_const, 0.0, 0.0);
  linear_const_penalty point_penalty = mvnormal_default_point_penalty(n, p, b_point);
  cor_mvnormal_capa_state S(n, p, std::move(precision_obj), penalty, point_penalty,
                            min_seg_len, max_seg_len);

  arma::mat x_t = x.t();  // Functions below expect column vectors.

  // run capa on all the data
  for (int i = 0; i < n; i++)
  {
    S = capa(
          capa_cor_mvnormal(
            update_cpts(
              prune(
                op(
                  cor_mvnormal_savings(
                    update_candidate_cpts(
                      update_vec_moments(
                        update_observation(std::move(S), x_t.col(i))
    ))))))));
  }


  // show the results
  auto locs_and_vars = mvcapa_locations_and_variables(S);
  // std::cout << "***************" << std::endl;
  // for(auto& lv : locs_and_vars)
  //   {
  //     std::cout << "(" << std::get<0>(lv)  << "," << std::get<1>(lv) << "," << std::get<2>(lv) << ")" << std::endl;
  //   }

  return format_capacc_output(locs_and_vars);
}
