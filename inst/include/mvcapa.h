
#ifndef ___MVCAPA_H___
#define ___MVCAPA_H___

#include "capa.h"
#include <vector>
#include <iostream>

namespace ostats
{

  struct mvcapa_state : public capa_state
  {
    std::vector<std::vector<int>> J;
    std::vector<std::vector<int>> J_collective;
    std::vector<std::vector<int>> J_point;
    mvcapa_state(const int&, const int&, const int&);
  };

  template <class state_type>
    std::list<std::tuple<int,int,int,double>> mvcapa_locations_and_variables(state_type S)
    {
      std::list<std::tuple<int,int,int,double>> lv;
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
      	  int cpt_dist = std::get<0>(S.anoms[pos]);
      	  int cpt = pos + 1 - cpt_dist;
      	  for (auto& j : S.J_collective[pos])
      	  {
      	    double change_size = (S.M1[pos + 1](j - 1) - S.M1[cpt](j - 1)) / cpt_dist;
    	      lv.push_back(std::make_tuple(cpt + 1, pos + 1, j, change_size));
      	  }
    	    pos = pos - cpt_dist;
    	  }
      	else
    	  {
      	  // Returns point-anomaly with indices starting from 1.
      	  for (auto& j : S.J_point[pos])
      	  {
      	    double change_size = S.M1[pos + 1](j - 1) - S.M1[pos](j - 1);
    	      lv.push_back(std::make_tuple(pos + 1, pos + 1, j, change_size));
      	  }
    	    pos -= 1;
    	  }
      }
    return lv;
    }

} // namespace ostats


#endif // ___MV_OP_H___
