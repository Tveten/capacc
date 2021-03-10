
#include "mvcapa.h"

namespace ostats
{
  mvcapa_state::mvcapa_state(const int& w, const int& min_sl, const int& max_sl)
  : capa_state(w, min_sl, max_sl), J(w)
  {
    J_collective.reserve(w);
    J_point.reserve(w);
  }
} // namespace ostats
