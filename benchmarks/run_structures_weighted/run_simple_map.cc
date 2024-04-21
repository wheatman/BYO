

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

#include "../run_weighted.h"
using weight_type = int;
#ifdef STD_MAP
#include <map>
static constexpr bool prefer_sorted_input = true;
using map_type = std::map<gbbs::uintE, weight_type>;
#endif

#ifdef STD_UNORDERED_MAP
#include <unordered_map>
static constexpr bool prefer_sorted_input = false;
using map_type = std::unordered_map<gbbs::uintE, weight_type>;
#endif

#ifdef ABSL_BTREE
#include "absl/container/btree_map.h"
static constexpr bool prefer_sorted_input = true;
using map_type = absl::btree_map<gbbs::uintE, weight_type>;
#endif

#ifdef ABSL_FLAT
#include "absl/container/flat_hash_map.h"
static constexpr bool prefer_sorted_input = false;
using map_type = absl::flat_hash_map<gbbs::uintE, weight_type>;
#endif

#ifdef USE_INPLACE
static constexpr bool use_inplace = true;
#else
static constexpr bool use_inplace = false;
#endif

using sym_graph_impl =
    gbbs::graph_implementations::symmetric_set_graph<map_type, weight_type,
                                                     /* inplace */ use_inplace,
                                                     prefer_sorted_input>;

using asym_graph_impl =
    gbbs::graph_implementations::asymmetric_set_graph<map_type, weight_type,
                                                      /* inplace */ use_inplace,
                                                      prefer_sorted_input>;

using graph_api = gbbs::full_api;

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
      using graph_t =
          gbbs::Graph<sym_graph_impl, /* symmetric */ true, graph_api>;
      auto G = gbbs::gbbs_io::read_weighted_symmetric_graph<graph_t>(
          iFile, mmap, binary);
      run_all(G, options);
    } else {
      using graph_t =
          gbbs::Graph<asym_graph_impl, /* symmetric */ false, graph_api>;
      auto G = gbbs::gbbs_io::read_weighted_symmetric_graph<graph_t>(
          iFile, mmap, binary);
      run_all(G, options);
      return -1;
    }
  }
  return 1;
}