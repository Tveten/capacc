

#include <iostream>


#include <fstream>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <set>

#include <limits>
#include <boost/circular_buffer.hpp>
#include <tuple>

#include <list>

typedef boost::circular_buffer<double> data_type;
typedef boost::circular_buffer<double> cost_type;
typedef boost::circular_buffer<double> penalty_type;
typedef boost::circular_buffer<int> cpts_type;
typedef std::tuple<int,int> anom_type;
typedef boost::circular_buffer<anom_type> anoms_type;



#include "median.h"
#include "mad.h"
#include "robustscale.h"
#include "observation.h"
#include "moments.h"
#include "op.h"
#include "changepoint.h"
#include "normal.changepoint.h"
#include "capa.h"
#include "normal.capa.h"
#include "normal.cost.h"
using namespace ostats;


#include <Eigen/Dense>
using Eigen::MatrixXd;


struct data_state
{
  boost::circular_buffer<double> data;
  data_state(const int& w);
};
data_state::data_state(const int& w) : data(w) {}


struct slope_state : data_state , changepoint_state
{
  slope_state(const int&,const double&);
}; 
slope_state::slope_state(const int& w, const double& penalty) : data_state(w) , changepoint_state(w,penalty) {};

template <class state_type>
state_type slope(state_type S)
{
  auto slope_cost = [&S](const int& i, const int& j)
    {
      if(j - i < 2)
	{
	  return std::numeric_limits<double>::infinity(); 
	} 
      MatrixXd X(j-i+1,2);
      MatrixXd y(j-i+1,1);
      int pos = 0;
      for(int k = i; k <= j; k++)
	{
	  X(pos,0) = double(pos);
	  X(pos,1) = 1.0;
	  y(pos,0) = S.data[k];
	  pos++;
	}
      auto Beta = (X.transpose()*X).inverse()*X.transpose()*y;
      auto y_hat = X*Beta;
      auto R = ((y - y_hat).transpose())*(y - y_hat);
      return R(0,0)/double(j-i+1);
      // return R(0,0);
    };
  
  auto n = S.data.size();
  S.C.clear();
  for(int i = 0; i < n; i++)
    {
      S.C.push_back(slope_cost(i,n-1));
    }
  
  return(S); 
}


template <class state_type>
state_type update_data(state_type S,const double& x)
{
  S.data.push_back(x);
  return(S);
}

int main()
{


  slope_state S(1000,3*log(200.0));

  
  
  std::ifstream ifile("./ramp-no-noise.csv", std::ios::in);
  
  //check to see that the file was opened correctly:
  if (!ifile.is_open()) {
    std::cerr << "There was a problem opening the input file!\n";
    exit(1);//exit or do additional error checking
  }

  std::vector<double> data;
  double num = 0.0;
  //keep storing values from the text file so long as data exists:
  while (ifile >> num) {
    data.push_back(num);
  }



  

  for(auto& x : data)
    {
      S = update_cpts(op(slope(update_penalty(update_data(std::move(S),x)))));
    }
  
  
  int count = 0;
  for(auto& cpt : S.cpts)
    {
      std::cout << count << " : " << cpt << std::endl;
      count++;
    }
    
  auto locs = changepoint_locations(S);

  std::cout << "***************" << std::endl;
  for(auto& loc : locs)
    {
      std::cout << loc << std::endl;
    }

  
  
  return 0;


  
}
