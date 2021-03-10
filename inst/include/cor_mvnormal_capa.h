
#ifndef ___COR_MVNORMAL_CAPA_H___
#define ___COR_MVNORMAL_CAPA_H___

#include "cor_mvnormal_saving.h"
#include "mvcapa.h"
#include "penalty.h"
#include "precision.h"

#include <cmath>

namespace ostats
{

  struct cor_mvnormal_capa_state : public mvcapa_state , cor_mvnormal_saving_state
  {
    linear_const_penalty point_penalty;
    cor_mvnormal_capa_state
    (
      const int&,
      const int&,
      precision&&,
      linear_const_penalty&,
      linear_const_penalty&,
      const int&,
      const int&
    );
  };

 template <class state_type>
   state_type capa_cor_mvnormal(state_type S)
   {
     auto n = S.F.size();
     // Collective anomaly
     S.C1 = S.F[n-1];
     S.J_collective.push_back(S.J[S.cpt]);

     // No anomaly
     S.C2 = S.F[n-2];

     // Point anomaly
     bqp::BQP_res point_saving = cor_mvnormal_saving(S.x, 1, S.prec_obj, S.point_penalty);
     S.C3 = S.F[n-2] - point_saving.max_value;
     S.J_point.push_back(point_saving.max_subset);
     return S;
   }

}  // namespace ostats

#endif // ___COR_MVNORMAL_CAPA_H___
