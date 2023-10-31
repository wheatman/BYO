

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


#define NO_TLX
#include "PMA/CPMA.hpp"

template <class T> class PMAWrapper {
  CPMA<PMA_SETTINGS<T>> pma;

public:
  PMAWrapper(auto begin, auto end) {
    for (auto it = begin; it != end; ++it) {
      pma.insert(*it);
    }
  }
  PMAWrapper() = default;
  template <class F> void map(F f) const {
    pma.template map<true>([&](auto el) { return f(el, {}); });
  }
  template <class F> void map_early_exit(F f) const {
    pma.template map<false>([&](auto el) { return f(el, {}); });
  }
  template <class F> void parallel_map(F f) const {
    pma.parallel_map([&](auto el) { return f(el, {}); });
  }
  void insert(auto el) { pma.insert(el); }
  void erase(auto el) { pma.remove(el); }
  void insert_sorted_batch(auto es, size_t n) { pma.insert_batch(es, n, true); }
  void remove_sorted_batch(auto es, size_t n) { pma.remove_batch(es, n, true); }

  size_t size() const { return pma.size(); }
};

#include "../run_unweighted.h"

using graph_impl = gbbs::graph_implementations::symmetric_set_graph<
    gbbs::symmetric_vertex, gbbs::empty, PMAWrapper<gbbs::uintE>>;

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
    std::cerr << "does not read in compressed files\n";
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