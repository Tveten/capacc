
#include "op.h"

namespace ostats
{
  // Changed: F(w) -> F(w + 1) because both start from index 0.
  op_state::op_state(const int& w, const int& min_sl, const int& max_sl)
  : F(w + 1) , T(w + 1), C(w) , min_seg_len(min_sl), max_seg_len(max_sl)
  {
    T.push_back(0);
    // F.push_back(-beta);
  }
} // namespace ostats
