
#include "capa.h"

namespace ostats
{
  capa_state::capa_state(const int& w, const double& penalty,const double& _beta_tilde) : beta_tilde(_beta_tilde) , anoms(w) , changepoint_state(w,penalty) {}
}
