
#include "changepoint.h"


namespace ostats
{
  changepoint_state::changepoint_state(const int& w,const double& _penalty) : op_state(w,_penalty) , cpts(w) , penalty(_penalty) {}
} // namespace ostats
