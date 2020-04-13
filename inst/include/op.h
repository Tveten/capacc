
#ifndef ___OP_H___
#define ___OP_H___

#include <boost/circular_buffer.hpp>
#include <vector>

namespace ostats
{
  struct op_state
  {
    boost::circular_buffer<double> F;
    boost::circular_buffer<int> T;
    // Changed so that C contains penalised costs. Needs to be this way to
    // incorporate multivariate costs and savings, because in this case
    // the penalised cost/savings are optimised as one entity. I.e., they
    // cannot be separated into distinct objects. Thus, I have moved the
    // penalty into the (penalised) cost objects.
    std::vector<double> C;
    // boost::circular_buffer<double> B;  Removed. See comment above.

    // Needed during pruning. And for multivariate cpt or anomaly detection,
    // it will not always be just the beta variable as in the univar case.
    // Initialised in the normal_capa_state or normal_changepoint_state
    // constructors, together with the first element of F.
    double pruning_const;

    int cpt;
    // Changed: Added min_distance as parameter of OP, rather than the cost,
    // because it is a parameter of OP that can be set independently of the cost.
    int min_seg_len;
    int max_seg_len;
    // TODO: Add max_distance
    op_state(const int&, const int&, const int&);
  };

  template <class state_type>
    state_type update_candidate_cpts(state_type S)
    {
      int n = S.F.size();
      if (n >= 2 * S.min_seg_len)
      {
        S.T.push_back(n - S.min_seg_len);
      }
      if (n - S.T.front() > S.max_seg_len)
      {
        S.T.pop_front();
      }
      // When the buffer has reached full circle, the indices must be moved
      // backwards accordingly.
      // if (S.T.back() == S.T.capacity() - S.min_distance + 1)
      // {
      //   for (auto& t : S.T)
      //   {
      //     t -= 1;
      //   }
      // }
      return std::move(S);
    }

  template <class state_type>
    state_type op(state_type S)
    {
      auto min_val = std::numeric_limits<double>::infinity();
      auto min_pos = 0;
      for(auto& t : S.T)
      {
        auto val = S.F[t] + S.C[t];
        if(val < min_val)
        {
          min_val = val;
          min_pos = t;
        }
      }
      S.F.push_back(min_val);
      S.cpt = min_pos;

      // auto Ffirst = S.F[0];
      // for(auto& f : S.F)
      // {
      //   f -= Ffirst;
      // }
      return S;
    }

  template <class state_type>
    state_type prune(state_type S)
    {
      for (auto it = S.T.begin(); it != S.T.end(); )
      {
        // std::cout << S.F.size() - 1 << " " << *it << ": " << S.F[*it] + S.C[*it] - S.beta << "  ------  " << S.F.back() << std::endl;
        // - S.pruning_const here because S.C is now the penalised cost.
    	  if(S.F[*it] + S.C[*it] - S.pruning_const > S.F.back())
    	  {
    	    // std::cout << S.F.size() - 1 << ", prune: " << *it << std::endl;
    	    it = S.T.erase(it);
    	  }
    	  else
    	  {
    	    it++;
    	  }
      }
      return std::move(S);
    }

} // namespace ostats



#endif // ___OP_H___
