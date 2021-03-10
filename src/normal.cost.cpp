
#include "normal.cost.h"


namespace ostats
{
  normal_cost_state::normal_cost_state
  (
    const constant_penalty& penalty_,
    const int& w
  ) :
    moments_state(w),
    penalty(penalty_) {}
} // namespace ostats
