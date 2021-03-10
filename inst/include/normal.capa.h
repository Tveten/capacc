
#ifndef ___NORMAL_CAPA_H___
#define ___NORMAL_CAPA_H___

#include "penalty.h"
#include "normal.cost.h"
#include "capa.h"

#include <cmath>

namespace ostats
{

  struct normal_capa_state : public capa_state , normal_cost_state
  {
    // double beta_tilde;
    constant_penalty point_penalty;
    normal_capa_state(const int&,
                      const constant_penalty&,
                      const constant_penalty&,
                      const int&,
                      const int&);
  };


 template <class state_type>
   state_type capa_normal_mean_var(state_type S)
   {
     auto n = S.F.size();
     S.C1 = S.F[n-1];
     S.C2 = S.F[n-2] + S.x*S.x;
     S.C3 = S.F[n-2] + 1 + log(S.x*S.x + 2*exp(-2*S.point_penalty.alpha_const)) + S.point_penalty.alpha_const;
     return(S);
   }

}


#endif // ___NORMAL_CAPA_H___
