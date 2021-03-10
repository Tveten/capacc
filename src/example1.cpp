

#include "ostats.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>


// form union of existing state models
struct sequential_mad_state : public mad_state , observation_state
{
  sequential_mad_state(const int&);
};

sequential_mad_state::sequential_mad_state(const int& w) : mad_state(w){}


// [[Rcpp::export]]
int main1()
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

    sequential_mad_state S(300);



  // run capa on all the data
  for(auto& x : data)
    {
      S = mad(median(update_observation(std::move(S),x)));
      std::cout << S.median << " : " << S.med_abs_dev << std::endl;
    }

  return 0;
}
