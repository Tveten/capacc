
#ifndef ___BQP_H___
#define ___BQP_H___

#include <RcppArmadillo.h>
#include <memory>
#include <list>
#include <iostream>
#include <algorithm>
#include <boost/functional/hash.hpp>


namespace bqp {
  struct node
  {
    std::weak_ptr<node> parent;
    std::shared_ptr<node> left_child;
    std::shared_ptr<node> right_child;
    char indicator;
    int level;
    int n_variables_included;
    double value;

    node(const std::weak_ptr<node>& parent_,
         const char& indicator_,
         const int& level_,
         const int& n_variables_included_,
         const double& value_);
  };

  typedef std::shared_ptr<node> node_ptr;
  typedef std::weak_ptr<node> weak_node_ptr;
  typedef std::vector<node_ptr> node_ptr_vec;
  typedef std::unordered_map<std::string, node_ptr_vec> umap_str_nodevec;
  typedef std::list<node_ptr> node_ptr_list;
  typedef std::vector<node_ptr> node_ptr_vec;
  typedef std::unordered_map<std::string, node_ptr_list> umap_str_nodelist;

  struct BQP {
    const arma::sp_mat A;
    const arma::mat b;
    const double c;
    const double k_star;
    const std::vector<std::vector<int>> nbs;
    const std::vector<std::vector<int>> extended_nbs;

    BQP(const arma::sp_mat&,
        const arma::mat&,
        const double&,
        const double&,
        const std::vector<std::vector<int>>&,
        const std::vector<std::vector<int>>&);
  };

  struct BQP_res {
    double max_value;
    std::vector<int> max_subset;
    int min_subset_size;
    BQP_res(const double& max_value_ = 0.0,
                     const std::vector<int>& max_subset_ = std::vector<int>(),
                     const int& min_subset_size_ = 0)
                   : max_value(max_value_),
                     max_subset(max_subset_),
                     min_subset_size(min_subset_size_) {}
  };

  BQP_res BQP_optim(const BQP&);

  std::vector<int> seq_int_reverse(const int& start, const int& end);

  arma::uvec get_position_vec(const std::vector<int>& levels,
                              const node_ptr& leaf_ptr);

  std::string get_position_at_levels(const std::vector<int>& levels,
                                     const node_ptr& leaf_ptr);

  void add_left_child(const node_ptr& parent);

  void add_right_child(const node_ptr& parent, const double& value);

  double right_child_penalised_savings(const node_ptr& parent,
                                       const int& d,
                                       const BQP& curr_BQP);

  node_ptr_vec grow_tree(const int& d,
                         const node_ptr_vec& parents,
                         const BQP& curr_BQP);


  node_ptr which_node_max(const node_ptr_vec& node_vec);

  node_ptr which_node_max(const node_ptr_list& node_list);

  node_ptr_vec select_parents(const int& d,
                              const node_ptr_vec& leaf_vec,
                              const std::vector<std::vector<int>>& extended_nbs);

  node_ptr init_tree();

  node_ptr_vec init_leaf_vec(const node_ptr& root);

  int min_variables_included(const node_ptr_vec& node_vec);

}  // namespace bqp

#endif // ___BQP_H___
