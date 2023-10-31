

// Usage:
// numactl -i all ./run_unweighted -src 10012 -s -m -rounds 3 twitter_SJ
// flags:
//   required:
//     -src: the source to compute the BFS from
//   optional:
//     -rounds : the number of times to run the algorithm
//     -c : indicate that the graph is compressed
//     -m : indicate that the graph should be mmap'd
//     -s : indicate that the graph is symmetric
//     -d : dump the output arrays to files, useful for debugging

#include "cpam/cpam.h"
#include "gbbs/bridge.h"

template <class T>
class CPAMWrapper {
  struct edge_entry {
    using key_t = T;  // a vertex_id
    static inline bool comp(key_t a, key_t b) { return a < b; }
  };
  #ifdef CPAM_COMPRESSED
  using edge_tree = cpam::diff_encoded_set<edge_entry, 64>;
  #else
  using edge_tree = cpam::pam_set<edge_entry, 64>;
  #endif
  edge_tree tree;

 public:
  size_t size() const { return tree.size(); }

  CPAMWrapper(auto begin, auto end) {
    for (auto it = begin; it != end; ++it) {
      tree.insert(*it);
    }
  }
  CPAMWrapper() = default;

  template <class F>
  void map(F f) const {
    auto map_f = [&](const auto &et) { f(et, {}); };
    ((edge_tree)tree).iterate_seq(map_f);
  }

  template <class F>
  void map_early_exit(F f) const {
    auto map_f = [&](const auto &et) { return !f(et, {}); };
    tree.foreach_cond(tree, map_f);
  }

  template <class F>
  void parallel_map(F f) const {
    auto map_f = [&](const auto &et, size_t i) { f(et, {}); };
    tree.foreach_index(tree, map_f);
  }

  // template <class F>
  // void parallel_map_early_exit(F f) const {
  // auto map_f = [&](const auto& et) { return f(et, {}); };
  // tree.foreach_cond_par(tree, map_f, []() { return true; });
  //}

  void insert(auto el) { tree.insert(el); }
  void erase(auto el) { tree.remove(el); }

  void insert_sorted_batch(auto es, size_t n) {
    auto replace = [](const auto &a, const auto &b) { return b; };
    tree = edge_tree::multi_insert_sorted(
        tree, std::ranges::subrange(es, es + n), replace);
  }
  void remove_sorted_batch(auto es, size_t n) {
    tree =
        edge_tree::multi_delete_sorted(tree, std::ranges::subrange(es, es + n));
  }
};

#include "../run_unweighted.h"

#ifdef USE_INPLACE
static constexpr bool use_inplace = true;
#else
static constexpr bool use_inplace = false;
#endif

using graph_impl = gbbs::graph_implementations::symmetric_set_graph<
    gbbs::symmetric_vertex, gbbs::empty, CPAMWrapper<gbbs::uintE>,
    /* inplace */ use_inplace>;

using graph_api = gbbs::full_api;

using graph_t = gbbs::Graph<graph_impl, /* symmetric */ true, graph_api>;

int main(int argc, char *argv[]) {
  gbbs::commandLine P(argc, argv, " [-s] <inFile>");
  char *iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");
  gbbs::run_all_options options;
  options.dump = P.getOptionValue("-d");
  options.rounds = P.getOptionLongValue("-rounds", 3);
  options.src = static_cast<gbbs::uintE>(P.getOptionLongValue("-src", 0));

  std::cout << "### Graph: " << iFile << std::endl;
  if (compressed) {
    std::cerr << "is always compressed, but reads in uncompressed files\n";
    return -1;
  } else {
    if (symmetric) {
      auto G = gbbs::gbbs_io::read_unweighted_symmetric_graph<graph_t>(
          iFile, mmap, binary);
      run_all<true>(G, options);
    } else {
      std::cerr << "does not support directed graphs yet\n";
      return -1;
    }
  }
  return 1;
}
