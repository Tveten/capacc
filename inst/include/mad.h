
#ifndef ___MAD_H___
#define ___MAD_H___

#include "median.h"
#include <cmath>

namespace ostats
{
  struct mad_state : median_state
  {
    double med_abs_dev;
    mad_state(const int&);
  };

  template <typename state_type>
    state_type mad(state_type S)
    {
      // 0 or 1 elements 
      if(S.X.size() < 2)
	{
	  S.med_abs_dev = 0.0;
	  return std::move(S);
	}
      auto it_lower = S.lower.begin();
      auto it_upper = S.upper.begin();
      int count = 0;
      auto n = S.lower.size() + S.upper.size();
      auto cut = n/2;
      do
	{
	  auto delta_lower = fabs(*it_lower - S.median);
	  auto delta_upper = fabs(*it_upper - S.median); 
	  if(delta_lower <= delta_upper)
	    {
	      S.med_abs_dev = delta_lower;
	      it_lower++;
	    }
	  else
	    {
	      S.med_abs_dev = delta_upper;
	      it_upper++;
	    }
	  count++;
	}
      while(count < cut);
      auto delta_lower = fabs(*it_lower - S.median);
      auto delta = fabs(*it_upper - S.median);
      if(delta_lower <= delta)
	{
	  delta = delta_lower;
	}
      // even numer of elements - interpolate
      if(n%2 == 0)
	{
	  S.med_abs_dev = (S.med_abs_dev + delta)/2.0;
	}
      else
	{
	  S.med_abs_dev = delta;
	}
      return std::move(S);
    }
  
}  // namespace ostats

# endif  // ___MAD_H___ 
