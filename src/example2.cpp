
#include "ostats.h"
#include <iostream>
#include <vector>
#include <limits>
#include <fstream>
#include <cmath>


  // struct normal_capa_mad_state : public normal_capa_state , mad_state
  // {
  //   normal_capa_mad_state(const int&, const double&, const double&, const int&);
  // };
  // normal_capa_mad_state::normal_capa_mad_state(const int& w, const double& beta, const double& beta_tilde, const int& m) : normal_capa_state(w,beta,beta_tilde,m) , mad_state(w) {}


// [[Rcpp::export]]
int main2()
{

  // read the data in
  std::ifstream ifile("./simple.test.dat", std::ios::in);
  std::vector<double> data;

  //check to see that the file was opened correctly:
  if (!ifile.is_open())
    {
      std::cerr << "There was a problem opening the input file!\n";
      exit(1);//exit or do additional error checking
    }

  double num;
  //keep storing values from the text file so long as data exists:
  while (ifile >> num)
    {
    data.push_back(num);
  }

  // add a point anomaly
    data[49] = 100.0;

  // initialise state
  double n = data.size();
  // double beta = 4.0*log(n);
  // double beta_tilde = 3.0*log(n);
  constant_penalty penalty(4.0 * log(n));
  constant_penalty point_penalty(3.0 * log(n));
  int min_seg_len = 2;
  int max_seg_len = std::numeric_limits<int>::infinity();
  normal_capa_state S(n, penalty, point_penalty, min_seg_len, max_seg_len);

  // run capa on all the data
  for(auto& x : data)
    {
      S = capa(
            capa_normal_mean_var(
              update_cpts(
                prune(
                  op(
                    normal_mean_var(
                      update_candidate_cpts(
                        update_moments(
                          update_observation(std::move(S), x)
    ))))))));
    }


  // show the results
  // auto locs = capa_locations(S);
  //
  // std::cout << "***************" << std::endl;
  // for(auto& loc : locs)
  //   {
  //     std::cout << "(" << std::get<0>(loc)  << "," << std::get<1>(loc) << ")" << std::endl;
  //   }

  return 0;
}
