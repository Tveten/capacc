
#ifndef ___NORMAL_COST_H___
#define ___NORMAL_COST_H___

#include "moments.h"
#include "penalty.h"
#include <cmath>
#include <limits>
#ifdef PARALLEL_COST_VECTOR
#include "omp.h"
#endif

namespace ostats
{
  struct normal_cost_state : moments_state
  {
    // Penalty moved to the cost structs. See comment in op.h.
    constant_penalty penalty;
    normal_cost_state(const constant_penalty&, const int&);
  };


  template <class state_type>
    state_type normal_mean(state_type S)
    {
      auto normal_mean_cost = [&S](const int& i, const int& j)
    	{
    	  auto dfirst = S.M1[j] - S.M1[i - 1];
    	  auto dsecond = S.M2[j] - S.M2[i - 1];
    	  return(dsecond-dfirst*dfirst/(j-i));
    	};
      auto n = S.M1.size() - 1;
      #ifdef _OPENMP
        #pragma omp parallel for
        for (int i = 0; i < S.T.size(); i++)
        {
          int t = S.T[i];
      	  S.C[t] = normal_mean_cost(t + 1, n) + S.penalty.alpha_const;  // Summing over obs (i + 1) to n.
        }
      #else
        for (auto& t : S.T) {
      	  S.C[t] = normal_mean_cost(t + 1, n) + S.penalty.alpha_const;  // Summing over obs (i + 1) to n.
        }
      #endif
      return std::move(S);
    }


  template <class state_type>
    state_type normal_mean_var(state_type S)
    {
      auto normal_mean_var_cost = [&S](const int& i, const int& j)
    	{
    	  auto m = j - i + 1;  // Number of observations to consider cost for.
    	  if(m < 2)
  	    {
    	    // Changed from returning Inf to Nan as Inf messes up the pruning.
    	    return std::numeric_limits<double>::quiet_NaN();
  	    }
    	  auto dfirst = S.M1[j] - S.M1[i - 1];
    	  auto dsecond = S.M2[j] - S.M2[i - 1];
    	  auto sigsq = ((dsecond)-(dfirst*dfirst)/m)/m;
    	  if(sigsq < 0.00000000000000001)
  	    {
  	      sigsq = 0.000000000000001;
  	    }
    	  return(m*(log(sigsq)+1));
    	};
      auto n = S.M1.size() - 1;
      #ifdef _OPENMP
        #pragma omp parallel for
        for (int i = 0; i < S.T.size(); i++)
        {
          int t = S.T[i];
      	  S.C[t] = normal_mean_var_cost(t + 1, n) + S.penalty.alpha_const;  // Summing over obs (i + 1) to n.
        }
      #else
        for (auto& t : S.T) {
      	  S.C[t] = normal_mean_var_cost(t + 1, n) + S.penalty.alpha_const;  // Summing over obs (i + 1) to n.
        }
      #endif
      return std::move(S);
    }

} // namespace ostats




#endif // ___NORMAL_COST_H___
