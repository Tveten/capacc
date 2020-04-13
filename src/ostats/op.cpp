
#include "op.h"

namespace ostats
{
  op_state::op_state(const int& w,const double& beta) : F(w) , C(w) , B(w) {F.push_back(-beta);}
} // namespace ostats
