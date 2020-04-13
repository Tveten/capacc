
#include "moments.h"

namespace ostats
{
  moments_state::moments_state(const int& w) : M1(w+1) , M2(w+1)
  {
    M1.push_back(0.0);
    M2.push_back(0.0);
  }
} // namespace ostats
