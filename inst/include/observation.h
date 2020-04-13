
#ifndef ___OBSERVATION_H___
#define ___OBSERVATION_H___

#include <RcppArmadillo.h>

namespace ostats
{
  struct observation_state
  {
    double x;
  };

  struct vec_observation_state
  {
    arma::vec x;
  };

  template <class state_type>
    state_type update_observation(state_type S, double x)
    {
      S.x = x;
      return S;
    }

  template <class state_type>
    state_type update_observation(state_type S, arma::vec x)
    {
      S.x = x;
      return S;
    }

} // namespace ostats


#endif // ___OBSERVATION_H___
