//
// Created by Yohsuke Murase on 2020/05/22.
//

#ifndef CPP__UNIONFIND_HPP_
#define CPP__UNIONFIND_HPP_

#include <cstdlib>
#include <vector>
#include <set>
#include <map>

class UnionFind {
 public:
  explicit UnionFind(size_t n) : parent(n) {
    for (size_t i = 0; i < n; i++) { parent[i] = i; }
  }
  size_t root(size_t i) {
    if (parent[i] != i) {
      size_t r = root(parent[i]);
      parent[i] = r;
    }
    return parent[i];
  }
  bool merge(size_t i, size_t j) {
    size_t ri = root(i);
    size_t rj = root(j);
    if (ri == rj) return false;
    else if (ri > rj) { parent[ri] = rj; }
    else if (ri < rj) { parent[rj] = ri; }
    return true;
  }
  std::map<size_t, std::set<size_t> > to_map() {
    std::map<size_t, std::set<size_t> > m;
    for (size_t i = 0; i < parent.size(); i++) {
      size_t r = root(i);
      m[r].insert(i);
    }
    return std::move(m);
  }
 private:
  std::vector<size_t> parent;
};

#endif //CPP__UNIONFIND_HPP_
