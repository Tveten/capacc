
#include "capa.h"

namespace ostats
{
  capa_state::capa_state(const int& w, const int& min_sl, const int& max_sl) : changepoint_state(w, min_sl, max_sl) , anoms(w) {}
}
