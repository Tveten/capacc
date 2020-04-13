
#include "normal.capa.h"

namespace ostats
{
  normal_capa_state::normal_capa_state(const int& w, const double& penalty,const double& _beta_tilde,const int& md) :  capa_state(w,penalty,_beta_tilde) , normal_cost_state(w,md) {}
} // namespace ostats
