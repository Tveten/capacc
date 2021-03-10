
#include "normal.changepoint.h"

namespace ostats
{
  normal_changepoint_state::normal_changepoint_state
  (
    const int& w,
    const double& beta,
    const int& min_sl,
    const int& max_sl
  ) :
    changepoint_state(w, min_sl, max_sl) , normal_cost_state(beta, w)
  {
    F.push_back(-beta);
    pruning_const = beta;
  }
} // namespace ostats
