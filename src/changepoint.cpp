
#include "changepoint.h"


namespace ostats
{
  // changepoint_state::changepoint_state(const int& w,const double& _beta, const int& md) : op_state(w, _beta, md) , cpts(w) , beta(_beta) {}
  changepoint_state::changepoint_state(const int& w, const int& min_sl, const int& max_sl)
  : op_state(w, min_sl, max_sl) , cpts(w) {}
} // namespace ostats
