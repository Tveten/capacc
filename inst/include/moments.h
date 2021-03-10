
#ifndef ___MOMENTS_H___
#define ___MOMENTS_H___

#include "observation.h"
#include <boost/circular_buffer.hpp>
#include <RcppArmadillo.h>
#include <iostream>

namespace ostats
{
  struct moments_state : observation_state
  {
    boost::circular_buffer<double> M1;
    boost::circular_buffer<double> M2;
    moments_state(const int&);
  };


  struct vec_moments_state : public vec_observation_state
  {
    boost::circular_buffer<arma::vec> M1;
    vec_moments_state(const int&, const int&);
  };

  template <class state_type>
    state_type update_moments(state_type S)
    {
      S.M1.push_back(S.M1.back() + S.x);
      S.M2.push_back(S.M2.back() + S.x*S.x);
    //   auto M1first = S.M1[0];
    //   for(auto& m1 : S.M1)
    // 	{
    // 	  m1 -= M1first;
    // 	}
    //   auto M2first = S.M2[0];
    //   for(auto& m2 : S.M2)
    // 	{
    // 	  m2 -= M2first;
    // 	}
      return S;
    }

  template <class state_type>
    state_type update_vec_moments(state_type S)
    {
      S.M1.push_back(S.M1.back() + S.x);

    //   auto M1first = S.M1[0];
    //   for(auto& m1 : S.M1)
    // 	{
    // 	  m1 -= M1first;
    // 	}
      return S;
    }

} // namespace ostats

#endif // ___MOMENTS_H___
