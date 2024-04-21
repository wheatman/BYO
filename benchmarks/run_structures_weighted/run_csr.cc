

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
    if (symmetric) {
      auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<int>(iFile, mmap);
      run_all(G, options);
    } else {
      auto G =
          gbbs::gbbs_io::read_compressed_asymmetric_graph<int>(iFile, mmap);
      run_all(G, options);
    }
  } else {
    if (symmetric) {
      using graph_impl =
          gbbs::graph_implementations::symmetric_graph<gbbs::symmetric_vertex,
                                                       int>;
      using graph_t = gbbs::Graph<graph_impl, /* symmetric */ true, graph_api>;
      auto G = gbbs::gbbs_io::read_weighted_symmetric_graph<graph_t>(
          iFile, mmap, binary);
      run_all(G, options);
    } else {
      using graph_impl =
          gbbs::graph_implementations::asymmetric_graph<gbbs::asymmetric_vertex,
                                                        int>;
      using graph_t = gbbs::Graph<graph_impl, /* symmetric */ false, graph_api>;
      auto G = gbbs::gbbs_io::read_weighted_asymmetric_graph<graph_t>(
          iFile, mmap, binary);
      run_all(G, options);
    }
  }
  return 1;
}