
#ifndef ___CHANGEPOINT_H___
#define ___CHANGEPOINT_H___

#include "op.h"

#include <boost/circular_buffer.hpp>
#include <list>

namespace ostats
{
  struct changepoint_state : public op_state
  {
    boost::circular_buffer<int> cpts;
    // double beta;
    // changepoint_state(const int&, const double&, const int&);
    changepoint_state(const int&, const int&, const int&);
  };

  template <class state_type>
    state_type update_cpts(state_type S)
    {
      // Not changepoints, but distance to change-point from current observation.
      int n = S.F.size() - 1;
      S.cpts.push_back(n - S.cpt);
      return S;
    }

  template <class state_type>
    std::list<int> changepoint_locations(state_type S)
    {
      std::list<int> loc;
      int pos = S.cpts.size() - 1;
      while(pos > 0)
    	{
        // Returns locations with indices starting from 1.
    	  loc.push_front(pos + 1 - S.cpts[pos]);
    	  pos = pos - S.cpts[pos];
    	}
      return loc;
    }
} // namespace ostats



#endif // ___CHANGEPOINT_H___
