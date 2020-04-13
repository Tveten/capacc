
#ifndef ___ROBUSTSCALE_H___
#define ___ROBUSTSCALE_H___

#include "median.h"
#include "mad.h"

namespace ostats
{
  template <typename state_type>
    state_type robustscale(state_type S)
    {
      S.x = (S.x - S.median)/(S.med_abs_dev*1.4826);
      return(S);
    }
} // namespace ostats




#endif // ___ROBUSTSCALE_H___
