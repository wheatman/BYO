

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

#include <functional>

#ifdef HOMEGROWN
#define PARLAY 1
#endif

#include "gbbs/bridge.h"
#include "terrace/terrace_graph.h"

#include "../run_unweighted.h"

template <class W> struct symmetric_terrace_graph {
  using weight_type = W;
  static constexpr bool binary = true;
  static_assert(binary);
  using vertex_weight_type = double;
  using edge_type = std::tuple<gbbs::uintE, W>;

  size_t M() const { return nodes.get_num_edges(); }

  size_t N() const { return nodes.get_num_vertices(); }

  auto degree(size_t i) const { return nodes.degree(i); }

  template <class F> void map_neighbors(size_t i, F f) const {
    if constexpr (binary) {
      W empty_weight = W();
      nodes.map_neighbors_no_early_exit(
          i, [&](auto src, auto dest) { return f(src, dest, empty_weight); });
    } else {
      nodes.map_neighbors_no_early_exit(i, f);
    }
  }
  template <class F> void map_neighbors_early_exit(size_t i, F f) const {
    if constexpr (binary) {
      W empty_weight = W();
      nodes.map_neighbors_early_exit(
          i, [&](auto src, auto dest) { return f(src, dest, empty_weight); });
    } else {
      nodes.map_neighbors_early_exit(i, f);
    }
  }
  template <class F> void parallel_map_neighbors(size_t i, F f) const {
    if constexpr (binary) {
      W empty_weight = W();
      nodes.parallel_map_neighbors_no_early_exit(
          i, [&](auto src, auto dest) { return f(src, dest, empty_weight); });
    } else {
      nodes.parallel_map_neighbors_no_early_exit(i, f);
    }
  }

  template <class F>
  void parallel_map_neighbors_early_exit(size_t i, F f) const {
    if constexpr (binary) {
      W empty_weight = W();
      nodes.parallel_map_neighbors_early_exit(
          i, [&](auto src, auto dest) { return f(src, dest, empty_weight); });
    } else {
      nodes.parallel_map_neighbors_early_exit(i, f);
    }
  }
  // ======================= Constructors and fields  ========================
  symmetric_terrace_graph() : vertex_weights(nullptr) {}
  // for now use the same constructor as the exsisting one for ease, but
  // probably write a more general constuctor for everyone to use later
  // TODO(wheatman) figure out what to do with _deletion_fn, probably should
  // just use unique pointer
  symmetric_terrace_graph(auto *v_data, size_t n, size_t m,
                          std::function<void()> _deletion_fn, edge_type *_e0,
                          vertex_weight_type *_vertex_weights = nullptr)
      : nodes(n), vertex_weights(_vertex_weights), deletion_fn(_deletion_fn) {
    printf("INITIALIZE TERRACE GRAPH WITH NODES %lu, EDGES %lu\n", n, m);
    printf("nodes in graph = %u\n", nodes.get_num_vertices());
    printf("edges in graph = %lu\n", nodes.get_num_edges());

    std::vector<uint32_t> srcs(m);
    std::vector<uint32_t> dests(m);

    printf("convert to edge list\n");
    gbbs::parallel_for(0, n, [&](uint32_t i) {
      for (size_t j = v_data[i].offset; j < v_data[i].offset + v_data[i].degree;
           j++) {
        srcs[j] = i;
        dests[j] = std::get<0>(_e0[j]);
        // TODO: weights if necessary
      }
    });

    printf("starting new build from batch\n");
    nodes.build_from_batch(srcs.data(), dests.data(), n, m);
    // nodes.add_edge_batch_no_perm(srcs.data(), dests.data(), m);

    printf("** STARTING VERIFY **\n");
    for (uint32_t i = 0; i < n; i++) {
      nodes.verify_neighbors(i, dests.data() + v_data[i].offset);
    }
    printf("** FINISHED VERIFY **\n");

    printf("AFTER INIT nodes in graph = %u\n", nodes.get_num_vertices());
    printf("AFTER INIT edges in graph = %lu\n", nodes.get_num_edges());
  }

  // TODO: write these
  // Move constructor
  symmetric_terrace_graph(symmetric_terrace_graph &&other) noexcept {
    printf("move constructor\n");
    /*
      nodes = other.nodes;
      vertex_weights = other.vertex_weights;
      other.vertex_weights = nullptr;
    */
  }
  /*
  // Move assignment
  symmetric_terrace_graph &
  operator=(symmetric_terrace_graph &&other) noexcept {
    nodes = other.nodes;
    vertex_weights = other.vertex_weights;
    other.vertex_weights = nullptr;
  }
  */
  // Copy constructor
  symmetric_terrace_graph(const symmetric_terrace_graph &other)
      : nodes(other.N()) {
    printf("copy constructor\n");
    for (size_t i = 0; i < N(); i++) {
      other.map_neighbors(
          i, [&](auto src, auto dest, auto val) { nodes.add_edge(src, dest); });
    }

    vertex_weights = nullptr;
    if (other.vertex_weights != nullptr) {
      vertex_weights = gbbs::new_array_no_init<vertex_weight_type>(N());
      gbbs::parallel_for(0, N(), [&](size_t i) {
        vertex_weights[i] = other.vertex_weights[i];
      });
    }
    deletion_fn = [&]() {
      if (vertex_weights != nullptr) {
        gbbs::free_array(vertex_weights, N());
      }
    };
  }
  ~symmetric_terrace_graph() { deletion_fn(); }

  void insert_sorted_batch(std::tuple<uint32_t, uint32_t> *es, size_t n) {
    nodes.add_edge_batch_no_perm(es, n);
  }

  void remove_sorted_batch(std::tuple<uint32_t, uint32_t> *es, size_t n) {
    parlay::parallel_for(0, n, [&](size_t i) {
      nodes.remove_edge(std::get<0>(es[i]), std::get<1>(es[i]));
    });
  }

  size_t get_memory_size() {
    return nodes.get_size();
  }

  // Graph Data
  graphstore::TerraceGraph nodes;
  vertex_weight_type *vertex_weights;
  // called to delete the graph
  std::function<void()> deletion_fn;
};

using graph_impl = symmetric_terrace_graph<gbbs::empty>;

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
  options.max_batch =
      static_cast<size_t>(P.getOptionLongValue("-max_batch", 1000000));
  options.src = static_cast<gbbs::uintE>(P.getOptionLongValue("-src", 0));
  options.inserts = P.getOptionValue("-i");

  std::cout << "### Graph: " << iFile << std::endl;
  if (compressed) {
    std::cerr << "does not support compression\n";
    return -1;
  } else {
    if (symmetric) {
      auto G = gbbs::gbbs_io::read_unweighted_symmetric_graph<graph_t>(
          iFile, mmap, binary);
      auto bytes_used = G.get_memory_size();
      std::cout << "total bytes used = " << bytes_used << "\n";
      run_all(G, options);
    } else {
      std::cerr << "does not support directed graphs yet\n";
      return -1;
    }
  }
  return 1;
}
