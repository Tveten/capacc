
#ifndef ___MEDIAN_H___
#define ___MEDIAN_H___

#include <set>
#include <boost/circular_buffer.hpp>
#include <limits>

namespace ostats
{

  struct median_state
  {
    std::multiset<double,std::greater<double> > lower;
    std::multiset<double,std::less<double> > upper;
    boost::circular_buffer<double> X;
    double median;
    median_state(const int&);
  };


  template <typename state_type>
    state_type median(state_type S)
    {
      auto x = S.x;
      // remove oldest element if full
      if(S.X.full())
	{
	  if(S.X[0] <= S.median)
	    {
	      auto it = S.lower.find(S.X[0]);
	      S.lower.erase(it); 
	    }
	  else
	    {
	      auto it = S.upper.find(S.X[0]);
	      S.upper.erase(it);
	    }
	  // rebalance if necessary
	  if(S.lower.size() ==  S.upper.size() + 2)
	    {
	      S.upper.insert(*(S.lower.begin()));
	      S.lower.erase(S.lower.begin());
	    }
	  else if(S.upper.size() ==  S.lower.size() + 2)
	    {
	      S.lower.insert(*(S.upper.begin()));
	      S.upper.erase(S.upper.begin());
	    }
	}
      S.X.push_back(x);
      // add new element
      auto n_lower = S.lower.size();
      auto n_upper = S.upper.size();
      // empty case
      if(n_lower == 0 && n_upper == 0)
	{
	  S.lower.insert(x);
	  S. median = x;
	  return std::move(S);
	}
      // heavy on low side
      if(n_lower == n_upper + 1)
	{
	  auto median = *(S.lower.begin());
	  if(x <= median)
	    {
	      S.upper.insert(median);
	      S.lower.erase(S.lower.begin());
	      S.lower.insert(x);
	    }
	  else
	    {
	      S.upper.insert(x);
	    }
	  median = *(S.upper.begin()) + *(S.lower.begin());
	  S.median = median/2.0;
	  return std::move(S);
	}
      // heavy on upper side
      if(n_upper == n_lower +1 )
	{
	  auto median = *(S.upper.begin());
	  if(x > median)
	    {
	      S.lower.insert(median);
	      S.upper.erase(S.upper.begin());
	      S.upper.insert(x);
	    }
	  else
	    {
	      S.lower.insert(x);
	    }
	  median = *(S.upper.begin()) + *(S.lower.begin());
	  S.median = median/2.0;
	  return std::move(S);
	}
      // non-empty and balanced
      auto median = *(S.upper.begin()) + *(S.lower.begin());
      median = median / 2.0;
      if(x <= median)
	{
	  S.lower.insert(x);
	  S.median = *(S.lower.begin());
	}
      else
	{
	  S.upper.insert(x);
	  S.median = *(S.upper.begin());
	}
      return std::move(S);
    }
  
  
} // namespace ostats


#endif  // ___MEDIAN_H___ 
