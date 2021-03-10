
#include "normal.cost.h"


namespace ostats
{
  normal_cost_state::normal_cost_state(const int& w, const int& _min_distance) : moments_state(w) , min_distance(_min_distance) {}
} // namespace ostats
