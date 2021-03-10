
#ifndef ___CAPA_H___
#define ___CAPA_H___

#include <tuple>
#include <boost/circular_buffer.hpp>
#include <list>
#include "changepoint.h"

namespace ostats
{
  struct capa_state : changepoint_state
  {
    boost::circular_buffer<std::tuple<int,int> > anoms;
    std::tuple<int,int> anom;
    double C1;
    double C2;
    double C3;
    capa_state(const int&, const int&, const int&);
  };


  template <class state_type>
    state_type capa(state_type S)
    {
      int n = S.cpts.size();
      if(S.C1 <= S.C2 && S.C1 <= S.C3)
    	{
    	  S.anoms.push_back(std::make_tuple(S.cpts[n - 1], 0));
    	  return S;
    	}
      if(S.C2 <= S.C3 && S.C2 <= S.C1)
    	{
    	  S.F[S.F.size() - 1] = S.C2;
    	  S.anoms.push_back(std::make_tuple(n - 1, 1));
    	  return S;
    	}
      S.F[S.F.size() - 1] = S.C3;
      S.anoms.push_back(std::make_tuple(n - 1, 2));
      return S;
    }

  template <class state_type>
    std::list<std::tuple<int,int> > capa_locations(state_type S)
    {
      std::list<std::tuple<int,int> > loc;
      int pos = S.anoms.size() - 1;
      while(pos > 0)
      {
    	  if(std::get<1>(S.anoms[pos]) == 1)
    	  {
    	    pos -= 1;
    	  }
      	else if(std::get<1>(S.anoms[pos]) == 0)
    	  {
      	  // Returns start- and end-points with indices starting from 1.
      	  int cpt = std::get<0>(S.anoms[pos]);
    	    loc.push_back(std::make_tuple(pos + 1 - cpt + 1, pos + 1));
    	    pos = pos - cpt;
    	  }
      	else
    	  {
      	  // Returns point-anomaly with indices starting from 1.
    	    loc.push_back(std::make_tuple(pos + 1,pos + 1));
    	    pos -= 1;
    	  }
      }
    return loc;
    }
} // namespace ostats



#endif // ___CAPA_H___
