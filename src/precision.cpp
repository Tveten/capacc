
#include "precision.h"

namespace ostats
{

  precision::precision
  (
    const arma::sp_mat& Q_,
    const std::vector<std::vector<int>>& nbs_,
    const std::vector<std::vector<int>>& extended_nbs_
  ) :
    Q(Q_),
    nbs(nbs_),
    extended_nbs(extended_nbs_) {}

  std::vector<std::vector<int>> lower_nbs(const arma::sp_mat& A)
  {
    std::vector<std::vector<int>> nbs;
    nbs.reserve(A.n_cols);
    std::vector<int> curr_nbs;
    nbs.push_back(curr_nbs);

    for (int j = 1; j < A.n_cols; j++)
    {
      std::vector<int> curr_nbs;
      for (int i = j - 1; i >= 0; i--)
      {
        if (A(i, j) != 0.0) curr_nbs.push_back(i + 1);
      }
      nbs.push_back(curr_nbs);
    }
    return nbs;
  }

  std::vector<std::vector<int>> extended_lower_nbs(std::vector<std::vector<int>> nbs)
  {
    // Assumes nbs are indexed starting from 1.

    int p = nbs.size();
    int band = 0;
    int curr_band;
    for (auto i = 1; i < p; i++)
    {
      if (nbs[i].size() > 0) curr_band = i + 1 - *min_element(nbs[i].begin(), nbs[i].end());
      else curr_band = 0;
      if (curr_band > band) band = curr_band;
    }

    for (auto i = 1; i < p - 1; i++)
    {
      for (auto j = i + 1; j < std::min(i + band, p); j++) {
        for (auto& k : nbs[j])
        {
          if (k < i + 1) nbs[i].push_back(k);
        }
      }
    }
    for (auto& curr_nbs : nbs)
    {
      std::sort(curr_nbs.begin(), curr_nbs.end(), std::greater<int>());
      curr_nbs.erase(unique(curr_nbs.begin(), curr_nbs.end()), curr_nbs.end());
    }
    return nbs;
  }

}


// [[Rcpp::export]]
std::vector<std::vector<int>> lower_nbs(const arma::sp_mat& A)
{
  return ostats::lower_nbs(A);
}

// [[Rcpp::export]]
std::vector<std::vector<int>> extended_lower_nbs(std::vector<std::vector<int>> nbs)
{
  return ostats::extended_lower_nbs(nbs);
}

// [[Rcpp::export]]
void test_precision_copy(arma::sp_mat& A)
{
  std::vector<std::vector<int>> nbs = lower_nbs(A);
  ostats::precision precision_obj(A, nbs, extended_lower_nbs(nbs));
  ostats::precision precision_obj2(precision_obj);
}

// [[Rcpp::export]]
void test_precision_move(arma::sp_mat& A)
{
  std::vector<std::vector<int>> nbs = lower_nbs(A);
  ostats::precision precision_obj(A, nbs, extended_lower_nbs(nbs));
  ostats::precision precision_obj2(std::move(precision_obj));
}
