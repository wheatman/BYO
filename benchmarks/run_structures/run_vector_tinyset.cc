

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


#ifdef HOMEGROWN
#define PARLAY 1
#endif


#include "SSTGraph/TinySet.hpp"
#include <limits>

class TinySetWrapper {
  SSTGraph::TinySetV<> ts;

public:
  TinySetWrapper() : ts(std::numeric_limits<uint32_t>::max()) {}
  TinySetWrapper(auto begin, auto end)
      : ts(std::numeric_limits<uint32_t>::max()) {
    for (auto it = begin; it != end; ++it) {
      ts.insert(*it);
    }
  }
  template <class F> void map(F f) const {
    ts.template map<true>([&](auto el) { return f(el, {}); });
  }
  template <class F> void map_early_exit(F f) const {
    ts.template map<false>([&](auto el) { return f(el, {}); });
  }

  template <class F> void parallel_map(F f) const {
    ts.template parallel_map<true>([&](auto el) { return f(el, {}); });
  }

  template <class F, typename Block_F = std::nullptr_t>
  void parallel_map_early_exit(F f, Block_F block_check = {}) const {
    ts.template parallel_map<false>([&](auto el) { return f(el, {}); });
  }

  void insert(auto el) { ts.insert(el); }
  void erase(auto el) { ts.remove(el); }

  size_t size() const { return ts.get_n(); }
};

#include "../run_unweighted.h"

using graph_impl = gbbs::graph_implementations::symmetric_set_graph<
    gbbs::symmetric_vertex, gbbs::empty, TinySetWrapper>;

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