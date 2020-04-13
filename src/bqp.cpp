
# include "bqp.h"

// namespace {
namespace bqp {
// using namespace bqp;

std::vector<int> seq_int_reverse(const int& start, const int& end)
{
  std::vector<int> v(end - start + 1);
  std::iota(v.rbegin(), v.rend(), start);
  return v;
}

std::vector<int> active_subset_at_levels(const std::vector<int>& levels,
                                         const node_ptr& leaf_ptr)
{
  // Returns position as indices of 1's only. Not 0-1's.
  // Assumes levels are sorted from largest to smallest.

  int n_levels = levels.size();
  std::vector<int> active_subset;
  active_subset.reserve(n_levels);
  node_ptr curr_node = leaf_ptr;
  int i = 0;
  while (i < n_levels)
  {
    if (curr_node->level == levels[i])
    {
      if (curr_node->indicator == '1')
      {
        active_subset.push_back(curr_node->level);
      }
      i += 1;
    }
    curr_node = curr_node->parent.lock();
  }
  return active_subset;
}


std::string position_at_levels(const std::vector<int>& levels,
                               const node_ptr& leaf_ptr)
{
  // Returns 0-1-string indicating position.
  int n_levels = levels.size();
  std::string position;
  position.reserve(n_levels);

  node_ptr curr_node = leaf_ptr;
  int i = 0;
  while (i < n_levels)
  {
    if (curr_node->level == levels[i])
    {
      position += curr_node->indicator;
      i += 1;
    }
    curr_node = curr_node->parent.lock();
  }
  return position;
}

void add_left_child(const node_ptr& parent)
{
  parent->left_child = node_ptr(std::make_shared<node>(
    parent, '0', parent->level + 1, parent->n_variables_included, parent->value
  ));
}

void add_right_child(const node_ptr& parent, const double& value)
{
  parent->right_child = node_ptr(std::make_shared<node>(
    parent, '1', parent->level + 1, parent->n_variables_included + 1, value
  ));
}

double right_child_penalised_savings(const node_ptr& parent,
                                     const int& d,
                                     const BQP& curr_BQP)
{
  std::vector<int> active_subset = active_subset_at_levels(curr_BQP.nbs[d - 1], parent);
  double sum_A_interactions = 0.0;
  for (auto& i : active_subset)
  {
    sum_A_interactions += curr_BQP.A(i - 1, d - 1);
  }
  return parent->value + curr_BQP.b(d - 1) + curr_BQP.A(d - 1, d - 1) + 2.0 * sum_A_interactions;
}

node_ptr_vec grow_tree(const int& d,
                       const node_ptr_vec& parents,
                       const BQP& curr_BQP)
{
  // Visitor that builds tree by splitting leaf nodes in node_list in a binary fashion.
  node_ptr_vec children;
  children.reserve(2 * parents.size());
    for(auto& parent : parents)
    {
      add_left_child(parent);
      add_right_child(parent, right_child_penalised_savings(parent, d, curr_BQP));
      children.push_back(parent->left_child);
      children.push_back(parent->right_child);
    }
    return children;
}


node_ptr which_node_max(const node_ptr_vec& node_vec)
{
  node_ptr max_node = node_vec.front();
  for (auto& curr_node : node_vec)
  {
    if (curr_node->value > max_node->value) max_node = curr_node;
  }
  return max_node;
}

node_ptr which_node_max(const node_ptr_list& node_list)
{
  node_ptr max_node = node_list.front();
  for(auto& curr_node : node_list)
  {
    if (curr_node->value > max_node->value) max_node = curr_node;
  }
  return max_node;
}

node_ptr_vec select_parents(const int& d,
                            const node_ptr_vec& leaf_vec,
                            const std::vector<std::vector<int>>& extended_nbs)
{
  // Visitor that selects which leaves in L that will become parent nodes.
  int curr_nbs_size = extended_nbs[d - 1].size();
  int prev_nbs_size = extended_nbs[d - 2].size();
  if (curr_nbs_size > prev_nbs_size) return leaf_vec;
  else
  {
    umap_str_nodevec grouped_leaves;
    grouped_leaves.reserve(std::pow(2, curr_nbs_size));
    for(auto& curr_leaf : leaf_vec)
    {
      std::string key = position_at_levels(extended_nbs[d - 1], curr_leaf);
      grouped_leaves[key].push_back(curr_leaf);
    }
    node_ptr_vec parents;
    parents.reserve(grouped_leaves.size());
    for (auto& leaf_group : grouped_leaves)
    {
      parents.push_back(which_node_max(leaf_group.second));
    }
    return parents;
  }
}

node_ptr init_tree()
{
  weak_node_ptr weak_null_ptr;
  node_ptr root(std::make_shared<node>(weak_null_ptr, ' ', 0, 0, 0));
  return root;
}

node_ptr_vec init_leaf_vec(const node_ptr& root)
{
  node_ptr_vec leaf_vec{root};
  return leaf_vec;
}

int min_variables_included(const node_ptr_vec& node_vec)
{
  int min_var = node_vec.front()->n_variables_included;
  for (auto& curr_node : node_vec)
  {
    if (curr_node->n_variables_included < min_var)
      min_var = curr_node->n_variables_included;
  }
  return min_var;
}

} // unnamed namespace



namespace bqp
{
node::node(const std::weak_ptr<node>& parent_,
           const char& indicator_,
           const int& level_,
           const int& n_variables_included_,
           const double& value_)
  : parent(parent_), indicator(indicator_), level(level_),
    n_variables_included(n_variables_included_), value(value_) {}

BQP::BQP(const arma::sp_mat& A_,
         const arma::mat& b_,
         const double& c_,
         const double& k_star_,
         const std::vector<std::vector<int>>& nbs_,
         const std::vector<std::vector<int>>& extended_nbs_)
  : A(A_), b(b_), c(c_), k_star(k_star_),
    nbs(nbs_), extended_nbs(extended_nbs_) {}

BQP_res BQP_optim(const BQP& curr_BQP)
{
  int p = curr_BQP.A.n_rows;
  node_ptr root = init_tree();
  node_ptr_vec leaf_vec = init_leaf_vec(root);
  leaf_vec = grow_tree(1, leaf_vec, curr_BQP);
  int min_subset_size= 0;
  int d = 2;
  while (d <= p && min_subset_size < curr_BQP.k_star)
  {
    leaf_vec = select_parents(d, leaf_vec, curr_BQP.extended_nbs);
    min_subset_size = min_variables_included(leaf_vec);
    leaf_vec = grow_tree(d, leaf_vec, curr_BQP);
    d++;
  }

  node_ptr max_node = which_node_max(leaf_vec);
  double max_value = max_node->value + curr_BQP.c;
  std::vector<int> max_subset = active_subset_at_levels(seq_int_reverse(1, d - 1), max_node);
  std::reverse(max_subset.begin(), max_subset.end());

  BQP_res results(max_value, max_subset, min_subset_size);
  return results;
}
} // namespace bqp

// void print_node_positions(const node_ptr_vec& node_vec)
// {
//   int curr_level = node_vec[0]->level;
//   std::vector<int> check_levels = seq_int_reverse(1, curr_level);
//
//   for(auto& curr_node : node_vec)
//   {
//     std::string position_str = position_at_levels(check_levels, curr_node);
//     std::cout << curr_node->value << " : " << position_str << std::endl;
//   }
// }
//
// void print_vector(const std::vector<int>& vec)
// {
//   for (int i = 0; i < vec.size(); i++)
//     std::cout << vec[i];
// }
//
// void print_uvector(const arma::uvec& vec)
// {
//   for (int i = 0; i < vec.n_elem; i++)
//     std::cout << vec(i);
// }
//
//   arma::uvec string_subset(const std::vector<int>& x, const std::string& binary_string)
//   {
//     arma::uvec subset_x = arma::zeros<arma::uvec>(x.size());
//     int i = 0;
//     for (auto& ind : binary_string)
//     {
//       if (ind == '1') subset_x(i) = x[i];
//       i++;
//     }
//     subset_x = subset_x.elem(find(subset_x > 0));
//     return subset_x;
//   }

// std::vector<int> seq_int(const int& start, const int& end)
// {
//   std::vector<int> vec(end - start + 1) ;
//   std::iota (std::begin(vec), std::end(vec), 1);
//   return vec;
// }
