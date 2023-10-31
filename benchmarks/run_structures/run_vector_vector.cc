

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

#include "gbbs/bridge.h"
#include <vector>

template <class T> struct vector_set {
  using weight_type = gbbs::empty;

  size_t size() const { return vec.size(); }

  template <class F> void map(F f) const {
    for (const auto &other : vec) {
      weight_type w{};
      f(other, w);
    }
  }

  template <class F> void map_early_exit(F f) const {
    for (const auto &other : vec) {
      weight_type w{};
      if (f(other, w)) {
        break;
      }
    }
  }

  template <class F> void parallel_map(F f) const {
    weight_type w{};
    gbbs::parallel_for(0, vec.size(), [&](auto i) { f(vec[i], w); });
  }

  template <class F, typename Block_F = std::nullptr_t>
  void parallel_map_early_exit(F f, Block_F block_check = {}) const {
    if (vec.size() < 1000) {
      map_early_exit(f);
    }
    size_t b_size = 2048;
    size_t n_blocks = vec.size() / b_size + 1;
    weight_type w{};
    gbbs::parallel_for(0, n_blocks, [&](size_t b) {
      if constexpr (!std::is_same_v<Block_F, std::nullptr_t>) {
        if (!block_check()) {
          return;
        }
      }
      size_t start = b * b_size;
      size_t end = std::min((b + 1) * b_size, static_cast<size_t>(vec.size()));
      for (size_t j = start; j < end; j++) {
        if (f(vec[j], w)) {
          return;
        }
      }
    });
  }

  vector_set() = default;
  vector_set(auto start, auto end) : vec(start, end) {}

private:
  gbbs::sequence<T> vec = {};
};

#include "../run_unweighted.h"

using graph_impl = gbbs::graph_implementations::symmetric_set_graph<
    gbbs::symmetric_vertex, gbbs::empty, vector_set<gbbs::uintE>>;

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
    std::cerr << "does not support compression\n";
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