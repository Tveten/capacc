
#ifndef ___NORMAL_CHANGEPOINT_H___
#define ___NORMAL_CHANGEPOINT_H___

#include "normal.cost.h"
#include "changepoint.h"
#include <limits>


namespace ostats
{

  struct normal_changepoint_state : public changepoint_state , normal_cost_state
    {
      normal_changepoint_state(const int&,const double&,const int&,const int&);
    };

} // namespace ostats


#endif // ___NORMAL_CHANGEPOINT_H___
