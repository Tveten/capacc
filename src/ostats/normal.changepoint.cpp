
#include "normal.changepoint.h"

namespace ostats
{
  normal_changepoint_state::normal_changepoint_state(const int& w, const double& penalty,const int& md) : changepoint_state(w,penalty) , normal_cost_state(w,md) {}
} // namespace ostats
