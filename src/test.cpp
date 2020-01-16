#include <stdio.h>
#include <iostream>
#include <memory>
#include <list>
#include <Rcpp.h>
using namespace Rcpp;

// void cast_test()
// {
//   int a = 129;
//
//   char b = a;
//
//   char c = (char)a;
//
//   int d = (int)b;
//
//   int e = a & 0xff;
//
//   // std::cout << a << " " << b << " " << c << " " << d << " " << e << std::endl;
// }

// void print_list(std::list<int> g)
// {
//   std::list<int>::iterator it;
//   for(it = g.begin(); it != g.end(); ++it)
//     std::cout << '\t' << *it;
// }

// [[Rcpp::export]]
void hash_map_test() {
  std::unordered_map<std::string, std::list<int>> umap;
  std::string key = "key";
  // std::list<int> int_list;
  // int_list.push_back(5);
  // umap[key] = int_list;
  umap[key].push_back(3);
    // grouped_leaves[key].push_back(node);

  std::cout << key << umap[key].front() << std::endl;
}
