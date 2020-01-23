#include <RcppArmadillo.h>
#include <memory>
#include <list>
#include <iostream>
#include <algorithm>
#include <boost/functional/hash.hpp>
using namespace Rcpp;

// NB - this requires c++14 or better. If your institutes computer systems do not have this intalled, tell your sytem
// administrator they are 5 years out of date !!!


namespace BQP_optim {
  struct node {
    std::shared_ptr<node> parent;
    std::shared_ptr<node> left_child;
    std::shared_ptr<node> right_child;
    char indicator;
    int level;
    int n_variables_included;
    double value;
  };

  typedef std::unordered_map<int, std::vector<int>> umap_int_vint;
  typedef std::shared_ptr<node> node_ptr;
  typedef std::list<node_ptr> node_ptr_list;
  typedef std::vector<node_ptr> node_ptr_vec;
  typedef std::unordered_map<std::string, node_ptr_list> umap_str_nodelist;
  typedef std::unordered_map<std::string, node_ptr_vec> umap_str_nodevec;

  struct BQP {
    arma::mat Q;
    arma::vec b;
    Rcpp::List penalty_components;
    umap_int_vint extended_nbs_map;

    BQP(const arma::mat& Q_,
        const arma::vec& b_,
        const Rcpp::List& penalty_components_,
        const Rcpp::List& extended_nbs_list_) {
      Q = Q_;
      b = b_;
      penalty_components = penalty_components_;
      int i = 0;
      for (std::vector<int> nbs : extended_nbs_list_) {
        extended_nbs_map[i] = nbs;
        i++;
      }
    }
  };
}; // namespace ubt

using namespace BQP_optim;

std::vector<int> seq_int_reverse(const int& start, const int& end) {
  int out_length = end - start + 1;
  std::vector<int> vec(out_length);
  for (int i = 0; i < out_length; i++) {
    vec[i] = end - i;
  }
  return vec;
}

arma::uvec get_position_vec(const std::vector<int>& levels,
                            const node_ptr& leaf_ptr) {
  int n_levels = levels.size();
  arma::uvec position = arma::zeros<arma::uvec>(n_levels);
  node_ptr curr_node = leaf_ptr;
  int i = 0;
  while (i < n_levels) {
    if (curr_node->level == levels[i]) {
      if (curr_node->indicator == '1') {
        position(i) = curr_node->level;
      }
      i += 1;
    }
    curr_node = curr_node->parent;
  }
  position = position.elem(find(position > 0));
  return position;
}

std::string get_position_at_levels(const std::vector<int>& levels,
                                   const node_ptr& leaf_ptr) {
  // Exception: levels must be smaller than or equal to leaf_ptr->level.
  // Exception: levels must be in decreasing order.
  int n_levels = levels.size();
  std::string position;
  position.reserve(n_levels);

  node_ptr curr_node = leaf_ptr;
  int i = 0;
  while (i < n_levels) {
    if (curr_node->level == levels[i]) {
      position += curr_node->indicator;
      i += 1;
    }
    curr_node = curr_node->parent;
  }
  return position;
}

void print_node_positions(const node_ptr_list& node_list) {
  int curr_level = node_list.front()->level;
  std::vector<int> check_levels = seq_int_reverse(1, curr_level);

  for(auto& curr_node : node_list) {
    std::string position_str = get_position_at_levels(check_levels, curr_node);
    std::cout << curr_node->value << " : " << position_str << std::endl;
  }
}

void print_vector(const std::vector<int>& vec) {
  for (int i = 0; i < vec.size(); i++)
    std::cout << vec[i];
}

void print_uvector(const arma::uvec& vec) {
  for (int i = 0; i < vec.n_elem; i++)
    std::cout << vec(i);
}


// [[Rcpp::export]]
arma::uvec string_subset(const std::vector<int>& x, const std::string& binary_string) {
  arma::uvec subset_x = arma::zeros<arma::uvec>(x.size());
  int i = 0;
  for (auto& ind : binary_string) {
    if (ind == '1') subset_x(i) = x[i];
    i++;
  }
  subset_x = subset_x.elem(find(subset_x > 0));
  return subset_x;
}

void add_left_child(const node_ptr& parent) {
  parent->left_child = node_ptr(new node);
  parent->left_child->parent = parent;
  parent->left_child->left_child = parent->left_child->right_child = nullptr;
  parent->left_child->indicator = '0';
  parent->left_child->n_variables_included = parent->n_variables_included;
  parent->left_child->level = parent->level + 1;
  parent->left_child->value = parent->value;
}

void add_right_child(const node_ptr& parent, const double& value) {
  parent->right_child = node_ptr(new node);
  parent->right_child->parent = parent;
  parent->right_child->left_child = parent->right_child->right_child = nullptr;
  parent->right_child->indicator = '1';
  parent->right_child->n_variables_included = parent->n_variables_included + 1;
  parent->right_child->level = parent->level + 1;
  parent->right_child->value = value;
}

// double right_child_penalised_savings(const node_ptr& parent, BQP& curr_BQP) {
double right_child_penalised_savings(const node_ptr& parent,
                                     const int& d,
                                     const std::vector<int>& nbs,
                                     const double& d_const_term,
                                     const arma::vec& Q_d,
                                     const double& beta) {
  // n must be multiplied with b and Q before sent here.
  double sum_Q_interactions = arma::accu(Q_d.elem(get_position_vec(nbs, parent) - 1));
  double additional_savings = d_const_term - 2.0 * sum_Q_interactions;
  return parent->value + additional_savings - beta;
}

node_ptr_list grow_tree(const int& d, const node_ptr_list& parents, BQP& curr_BQP) {
  // Visitor that builds tree by splitting leaf nodes in node_list in a binary fashion.
  arma::vec Q_d = curr_BQP.Q.col(d - 1);
  double d_const_term = 2.0 * curr_BQP.b(d - 1) - Q_d(d - 1);
  std::vector<int> nbs = curr_BQP.extended_nbs_map[d - 1];
  double beta = curr_BQP.penalty_components["beta"];
  node_ptr_list children;
  for(auto& parent : parents) {
    add_left_child(parent);
    add_right_child(parent, right_child_penalised_savings(parent, d, nbs, d_const_term, Q_d, beta));
    children.push_back(parent->left_child);
    children.push_back(parent->right_child);
  }
  return children;
}

node_ptr which_node_max(const node_ptr_list& node_list) {
  node_ptr max_node = node_list.front();
  for(auto& curr_node : node_list) {
    if (curr_node->value > max_node->value) max_node = curr_node;
  }
  return max_node;
}

node_ptr_list select_parents(const int& d, const node_ptr_list& leaf_list,
                             umap_int_vint& extended_nbs_map) {
  // Visitor that selects which leaves in L that will become parent nodes.
  node_ptr_list parents;
  if (extended_nbs_map[d - 1].size() > extended_nbs_map[d - 2].size()) {
    parents = leaf_list;
  }
  else {
    umap_str_nodelist grouped_leaves;
    for(auto& curr_leaf : leaf_list) {
      std::string key = get_position_at_levels(extended_nbs_map[d - 1], curr_leaf);
      grouped_leaves[key].push_back(curr_leaf);
    }
    for (auto& leaf_group : grouped_leaves) {
      parents.push_back(which_node_max(leaf_group.second));
    }
  }
  return parents;
}

node_ptr init_tree() {
  node_ptr root(new node);
  root->parent = root->left_child = root->right_child = nullptr;
  root->level = 0;
  root->n_variables_included = 0;
  root->value = 0;
  return root;
}

node_ptr_list init_leaf_list(const node_ptr& root) {
  node_ptr_list leaf_list;
  leaf_list.push_back(root);
  return leaf_list;
}

int min_variables_included(const node_ptr_list& node_list) {
  int min_var = node_list.front()->n_variables_included;
  for(auto& curr_node : node_list) {
    if (curr_node->n_variables_included < min_var)
      min_var = curr_node->n_variables_included;
  }
  return min_var;
}

// [[Rcpp::export]]
Rcpp::List optimise_savings(const arma::mat& Q,
                            const arma::vec& b,
                            const Rcpp::List& penalty_components,
                            const Rcpp::List& extended_nbs_list) {

  // Later: Create M from Q.
  int p = Q.n_rows;
  BQP curr_BQP(Q, b, penalty_components, extended_nbs_list);
  node_ptr root = init_tree();
  node_ptr_list leaf_list = init_leaf_list(root);
  leaf_list = grow_tree(1, leaf_list, curr_BQP);
  int k_star = penalty_components["k_star"];
  int min_n_variables_included = 0;
  int d = 2;
  while (d <= p && min_n_variables_included < k_star) {
    leaf_list = select_parents(d, leaf_list, curr_BQP.extended_nbs_map);
    min_n_variables_included = min_variables_included(leaf_list);
    leaf_list = grow_tree(d, leaf_list, curr_BQP);
    d++;
  }

  node_ptr max_node = which_node_max(leaf_list);
  std::vector<int> all_levels = seq_int_reverse(1, d - 1);
  arma::uvec J_max_uvec = get_position_vec(all_levels, max_node);
  std::vector<int> J_max = arma::conv_to<std::vector<int>>::from(J_max_uvec);
  double alpha = penalty_components["alpha"];
  double B_max = max_node->value - alpha;

  Rcpp::List results_list = List::create(_["B_max"] = B_max,
                                         _["J_max"] = J_max,
                                         _["k_min"] = min_n_variables_included);
  return results_list;
}

////////////// Vector ----------------------------------------------------------
void print_node_positions(const node_ptr_vec& node_vec) {
  int curr_level = node_vec[0]->level;
  std::vector<int> check_levels = seq_int_reverse(1, curr_level);

  for(auto& curr_node : node_vec) {
    std::string position_str = get_position_at_levels(check_levels, curr_node);
    std::cout << curr_node->value << " : " << position_str << std::endl;
  }
}

node_ptr_vec grow_tree(const int& d, const node_ptr_vec& parents, BQP& curr_BQP) {
  // Visitor that builds tree by splitting leaf nodes in node_list in a binary fashion.
  arma::vec Q_d = curr_BQP.Q.col(d - 1);
  double d_const_term = 2.0 * curr_BQP.b(d - 1) - Q_d(d - 1);
  std::vector<int> nbs = curr_BQP.extended_nbs_map[d - 1];
  double beta = curr_BQP.penalty_components["beta"];
  node_ptr_vec children(2 * parents.size(), nullptr);
  int i = 0;
  for(auto& parent : parents) {
    add_left_child(parent);
    add_right_child(parent, right_child_penalised_savings(parent, d, nbs, d_const_term, Q_d, beta));
    children[i] = parent->left_child;
    children[i + 1] = parent->right_child;
    i += 2;
  }
  return children;
}

node_ptr which_node_max(const node_ptr_vec& node_vec) {
  node_ptr max_node = node_vec[0];
  for (int i = 1; i < node_vec.size(); i ++) {
    if (node_vec[i]->value > max_node->value) max_node = node_vec[i];
  }
  return max_node;
}

node_ptr_vec select_parents(const int& d, const node_ptr_vec& leaf_vec,
                             umap_int_vint& extended_nbs_map) {
  // Visitor that selects which leaves in L that will become parent nodes.
  std::vector<int> curr_nbs = extended_nbs_map[d - 1];
  std::vector<int> prev_nbs = extended_nbs_map[d - 2];
  if (curr_nbs.size() > prev_nbs.size()) {
    return leaf_vec;
  }
  else {
    umap_str_nodelist grouped_leaves;
    for(auto& curr_leaf : leaf_vec) {
      std::string key = get_position_at_levels(curr_nbs, curr_leaf);
      grouped_leaves[key].push_back(curr_leaf);
    }
    node_ptr_vec parents(grouped_leaves.size(), nullptr);
    int i = 0;
    for (auto& leaf_group : grouped_leaves) {
      parents[i] = which_node_max(leaf_group.second);
      i++;
    }
    return parents;
  }
}

node_ptr_vec init_leaf_vec(const node_ptr& root) {
  node_ptr_vec leaf_vec{root};
  return leaf_vec;
}

int min_variables_included(const node_ptr_vec& node_vec) {
  int min_var = node_vec[0]->n_variables_included;
  for (int i = 1; i < node_vec.size(); i ++) {
    if (node_vec[i]->n_variables_included < min_var)
      min_var = node_vec[i]->n_variables_included;
  }
  return min_var;
}

// [[Rcpp::export]]
Rcpp::List optimise_savings_vec(const arma::mat& Q,
                                const arma::vec& b,
                                const Rcpp::List& penalty_components,
                                const Rcpp::List& extended_nbs_list) {

  // Later: Create M from Q.
  int p = Q.n_rows;
  BQP curr_BQP(Q, b, penalty_components, extended_nbs_list);
  int k_star = penalty_components["k_star"];
  node_ptr root = init_tree();
  node_ptr_vec leaf_vec = init_leaf_vec(root);
  leaf_vec = grow_tree(1, leaf_vec, curr_BQP);
  int min_n_variables_included = 0;
  int d = 2;
  while (d <= p && min_n_variables_included < k_star) {
    leaf_vec = select_parents(d, leaf_vec, curr_BQP.extended_nbs_map);
    min_n_variables_included = min_variables_included(leaf_vec);
    leaf_vec = grow_tree(d, leaf_vec, curr_BQP);
    d++;
  }

  node_ptr max_node = which_node_max(leaf_vec);
  std::vector<int> all_levels = seq_int_reverse(1, d - 1);
  arma::uvec J_max_uvec = get_position_vec(all_levels, max_node);
  std::vector<int> J_max = arma::conv_to<std::vector<int>>::from(J_max_uvec);
  double alpha = penalty_components["alpha"];
  double B_max = max_node->value - alpha;

  Rcpp::List results_list = List::create(_["B_max"] = B_max,
                                         _["J_max"] = J_max,
                                         _["k_min"] = min_n_variables_included);
  return results_list;
}

// Old -------------------------------------------------------------------------
void add_left_child_test(const std::shared_ptr<node>& parent, const double& val) {
  parent->left_child = std::shared_ptr<node>(new node);
  parent->left_child->parent = parent;
  parent->left_child->left_child = parent->left_child->right_child = nullptr;
  parent->left_child->indicator = '0';
  parent->left_child->level = parent->level + 1;
  parent->left_child->value = val;
}

void add_right_child_test(const std::shared_ptr<node>& parent, const double& val) {
  parent->right_child = std::shared_ptr<node>(new node);
  parent->right_child->parent = parent;
  parent->right_child->left_child = parent->right_child->right_child = nullptr;
  parent->right_child->indicator = '1';
  parent->right_child->level = parent->level + 1;
  parent->right_child->value = val;
}

std::vector<char> get_position(const std::shared_ptr<node>& leaf_ptr) {
  int curr_level = leaf_ptr->level;
  std::vector<char> position(curr_level, 0);
  std::shared_ptr<node> curr_node = leaf_ptr;
  for(int i = 0; i < curr_level; i++) {
    position[i] = curr_node->indicator;
    curr_node = curr_node->parent;
  }
  return position;
}

std::vector<char> get_position_after_level(const int start_level,
                                           const std::shared_ptr<node>& leaf_ptr) {
  // Exception: lead_ptr->level has to be greater than or equal to start_level.
  int curr_level = leaf_ptr->level;
  std::vector<char> position(curr_level - start_level + 1, 0);
  std::shared_ptr<node> curr_node = leaf_ptr;
  for(int i = 0; i < position.size(); i++) {
    position[i] = curr_node->indicator;
    curr_node = curr_node->parent;
  }
  return position;
}

// [[Rcpp::export]]
void list_test(Rcpp::List L) {
  for (auto& vec : L) {
    print_vector(vec);
  }
}

// [[Rcpp::export]]
std::vector<int> seq_int(const int& start, const int& end) {
  std::vector<int> v(end - start + 1);
  std::iota (std::begin(v), std::end(v), start);
  return v;
}

// [[Rcpp::export]]
std::unordered_map<std::string, double> test_map() {
  std::unordered_map<std::string, double> test_umap;
  test_umap["hehe"] = 2.0;
  test_umap["h"] = 1.0;
  return test_umap;
}

std::list<std::shared_ptr<node>> grow_tree_test(const std::list<std::shared_ptr<node>>& L) {
  // Visitor that builds tree by splitting leaf nodes in L in a binary fashion.
  std::list<std::shared_ptr<node>> new_L;
  double val = 0.0;
  for(auto& root : L) {
    add_left_child_test(root, val);
    val += 1.0;
    add_right_child_test(root, val);
    val += 1.0;
    new_L.push_back(root->left_child);
    new_L.push_back(root->right_child);
  }
  return new_L;
}

// [[Rcpp::export]]
void test_ubt(const int& n_levels = 4) {
  // create the root node
  std::shared_ptr<node> root(new node);
  // initialise root
  root->parent = root->left_child = root->right_child = nullptr;
  root->level = 0;
  // initialise visitor list
  std::list<std::shared_ptr<node>> L;
  L.push_back(root);
  L = grow_tree_test(L);
  L = grow_tree_test(L);
  print_node_positions(L);

  for (int d = 3; d <= n_levels; d++) {
    std::cout << "d = " << d << ":" << std::endl;
    std::vector<int> M_d = {d - 1, d - 2};
    // L = select_parents(L, M_d);
    print_node_positions(L);
    L = grow_tree_test(L);
    print_node_positions(L);
  }
}

// [[Rcpp::export]]
void which_max_test() {
  std::shared_ptr<node> root(new node);
  // initialise root
  root->parent = root->left_child = root->right_child = nullptr;
  root->level = 0;
  // initialise visitor list
  std::list<std::shared_ptr<node>> L;
  L.push_back(root);
  L = grow_tree_test(L);
  L = grow_tree_test(L);
  L = grow_tree_test(L);
  for(auto& root : L) {
    std::cout << root->value;
  }
  std::cout << std::endl;

  std::cout << which_node_max(L)->value << std::endl;
}

// [[Rcpp::export]]
void test_rev_vec() {
  std::vector<int> vec = seq_int_reverse(2, 5);
  print_vector(vec);
}

// [[Rcpp::export]]
double test_accu() {
  arma::uvec v{1, 2, 3};
  v = v.elem(find(v < 1));
  return arma::accu(v);
}

// [[Rcpp::export]]
Rcpp::List optimise_savings_vec_old(const arma::mat& Q,
                                    const arma::vec& b,
                                    const Rcpp::List& penalty_components,
                                    const Rcpp::List& extended_nbs_list) {

  // Later: Create M from Q.
  int p = Q.n_rows;
  BQP curr_BQP(Q, b, penalty_components, extended_nbs_list);
  node_ptr root = init_tree();
  node_ptr_vec leaf_vec = init_leaf_vec(root);
  leaf_vec = grow_tree(1, leaf_vec, curr_BQP);

  for (int d = 2; d <= p; d++) {
    leaf_vec = select_parents(d, leaf_vec, curr_BQP.extended_nbs_map);
    leaf_vec = grow_tree(d, leaf_vec, curr_BQP);
  }
  node_ptr max_node = which_node_max(leaf_vec);
  std::vector<int> all_levels = seq_int_reverse(1, p);
  arma::uvec J_max_uvec = get_position_vec(all_levels, max_node);
  std::vector<int> J_max = arma::conv_to<std::vector<int>>::from(J_max_uvec);
  double alpha = penalty_components["alpha"];
  double B_max = max_node->value - alpha;

  Rcpp::List results_list = List::create(_["B_max"] = B_max, _["J_max"] = J_max);
  return results_list;
}

// [[Rcpp::export]]
void test_grow_tree(const arma::mat& Q,
                    const arma::vec& b,
                    const Rcpp::List& penalty_components,
                    const Rcpp::List& extended_nbs_list) {
  BQP curr_BQP(Q, b, penalty_components, extended_nbs_list);
  node_ptr root = init_tree();
  node_ptr_list leaf_list = init_leaf_list(root);
  leaf_list = grow_tree(1, leaf_list, curr_BQP);
  leaf_list = grow_tree(2, leaf_list, curr_BQP);
  leaf_list = grow_tree(3, leaf_list, curr_BQP);
  leaf_list = grow_tree(4, leaf_list, curr_BQP);
  leaf_list = grow_tree(5, leaf_list, curr_BQP);
  leaf_list = grow_tree(6, leaf_list, curr_BQP);
  leaf_list = select_parents(7, leaf_list, curr_BQP.extended_nbs_map);
  leaf_list = grow_tree(7, leaf_list, curr_BQP);

  // print_node_positions(leaf_list);
}

// [[Rcpp::export]]
Rcpp::List optimise_savings_old(const arma::mat& Q,
                                const arma::vec& b,
                                const Rcpp::List& penalty_components,
                                const Rcpp::List& extended_nbs_list) {

  // Later: Create M from Q.
  int p = Q.n_rows;
  BQP curr_BQP(Q, b, penalty_components, extended_nbs_list);
  node_ptr root = init_tree();
  node_ptr_list leaf_list = init_leaf_list(root);
  leaf_list = grow_tree(1, leaf_list, curr_BQP);
  // print_node_positions(leaf_list);

  for (int d = 2; d <= p; d++) {
    leaf_list = select_parents(d, leaf_list, curr_BQP.extended_nbs_map);
    // print_node_positions(leaf_list);
    leaf_list = grow_tree(d, leaf_list, curr_BQP);
    // print_node_positions(leaf_list);
  }
  node_ptr max_node = which_node_max(leaf_list);
  std::vector<int> all_levels = seq_int_reverse(1, p);
  arma::uvec J_max_uvec = get_position_vec(all_levels, max_node);
  std::vector<int> J_max = arma::conv_to<std::vector<int>>::from(J_max_uvec);
  double alpha = penalty_components["alpha"];
  double B_max = max_node->value - alpha;

  int min_n_variables_included = 1;
  Rcpp::List results_list = List::create(_["B_max"] = B_max,
                                         _["J_max"] = J_max,
                                         _["k_min"] = min_n_variables_included);
  return results_list;
}
