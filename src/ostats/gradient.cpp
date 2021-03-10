
#include <iostream>
#include <Eigen/Dense>
using Eigen::MatrixXd;

#include <boost/circular_buffer.hpp>


int main()
{

  boost::circular_buffer<double> Y(1000);
  for(double y = 3.0; y < 20; y += 0.2)
    {
      Y.push_back(y);
    }
  boost::circular_buffer<double> T(1000);
  double t = 0.0;
  for(auto& y : Y)
    {
      T.push_back(t);
      t += 0.1;
    }

  MatrixXd X(Y.size(),2);
  int i = 0;
  for(auto& t : T)
    {
      X(i,0) = t;
      X(i,1) = 1.0;
      i++;
    }

  MatrixXd y(Y.size(),1);
  i = 0;
  for(auto& val : Y)
    {
      y(i,0) = val;
      i++;
    }
      
  auto Beta = (X.transpose()*X).inverse()*X.transpose()*y;
  auto y_hat = X*Beta;
  auto R = ((y - y_hat).transpose())*(y - y_hat);
  

  std::cout << R(0,0) << std::endl;
  
  
  MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  std::cout << m << std::endl;
}
