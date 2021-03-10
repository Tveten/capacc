
#include "normal.capa.h"

namespace ostats
{
  // normal_capa_state::normal_capa_state(const int& w, const double& beta,const double& _beta_tilde, const int& md) : capa_state(w, beta, beta_tilde, md) , normal_cost_state(w) {}
  normal_capa_state::normal_capa_state
  (
    const int& w,
    const constant_penalty& penalty_,
    const constant_penalty& point_penalty_,
    const int& min_sl,
    const int& max_sl
  ) :
    capa_state(w, min_sl, max_sl) , normal_cost_state(penalty_, w), point_penalty(point_penalty_)
  {
    F.push_back(-penalty.alpha_const);
    pruning_const = penalty.alpha_const;
  }
} // namespace ostats
