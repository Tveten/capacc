
#include "cor_mvnormal_capa.h"

namespace ostats
{
  cor_mvnormal_capa_state::cor_mvnormal_capa_state
  (
    const int& w,
    const int& p,
    precision&& precision_,
    linear_const_penalty& penalty_,
    linear_const_penalty& point_penalty_,
    const int& min_sl,
    const int& max_sl
  ) :
    mvcapa_state(w, min_sl, max_sl),
    cor_mvnormal_saving_state(w, p, std::move(precision_), penalty_),
    point_penalty(point_penalty_)
  {
    F.push_back(0.0);
    pruning_const = penalty_.P(precision_.Q.n_cols);
  }
} // namespace ostats
