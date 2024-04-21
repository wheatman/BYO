

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
#include <memory>
#include <random>

#include "../run_unweighted.h"
#include "gbbs/bridge.h"

#ifdef CSR_SHUFFLE
static constexpr bool csr_shuffle = true;
#else
static constexpr bool csr_shuffle = false;
#endif

template <class node_t, class edge_t, bool shuffle = false> struct CSR {
  using weight_type = gbbs::empty;
  using vertex_weight_type = double;
  size_t num_vertices() const { return n; }
  size_t num_edges() const { return m; }

  size_t degree(size_t i) const {
    return vertex_offsets[i + 1] - vertex_offsets[i];
  }

  template <class F> void map_neighbors(size_t i, F f) const {
    weight_type w{};
    for (edge_t j = vertex_offsets[i]; j < vertex_offsets[i + 1]; j++) {
      f(i, edges[j], w);
    }
  }

  template <class F> void map_neighbors_early_exit(size_t i, F f) const {
    weight_type w{};
    for (edge_t j = vertex_offsets[i]; j < vertex_offsets[i + 1]; j++) {
      if (f(i, edges[j], w)) {
        break;
      }
    }
  }

  template <class F> void parallel_map_neighbors(size_t i, F f) const {
    weight_type w{};
    gbbs::parallel_for(vertex_offsets[i], vertex_offsets[i + 1],
                       [&](auto j) { f(i, edges[j], w); });
  }

  template <class F, typename Block_F = std::nullptr_t>
  void parallel_map_neighbors_early_exit(size_t i, F f,
                                         Block_F block_check = {}) const {
    size_t size = degree(i);
    if (size < 1000) {
      map_early_exit(f);
    }
    size_t b_size = 2048;
    size_t n_blocks = size / b_size + 1;
    weight_type w{};
    gbbs::parallel_for(0, n_blocks, [&](size_t b) {
      if constexpr (!std::is_same_v<Block_F, std::nullptr_t>) {
        if (!block_check()) {
          return;
        }
      }
      size_t start = b * b_size;
      size_t end = std::min((b + 1) * b_size, static_cast<size_t>(size));
      for (size_t j = start; j < end; j++) {
        if (f(i, edges[j + vertex_offsets[i]], w)) {
          return;
        }
      }
    });
  }

  CSR() = default;
  CSR(auto *v_data, size_t n, size_t m, std::function<void()> _deletion_fn,
      auto *_e0, vertex_weight_type *_vertex_weights = nullptr)
      : n(n), m(m), vertex_weights(_vertex_weights), deletion_fn(_deletion_fn) {
    vertex_offsets.reset((edge_t *)malloc((n + 1) * sizeof(edge_t)));
    edges.reset((node_t *)malloc((m) * sizeof(node_t)));
    gbbs::parallel_for(0, n, [&](size_t i) {
      vertex_offsets[i] = v_data[i].offset;
      for (size_t j = 0; j < v_data[i].degree; j++) {
        edges[vertex_offsets[i] + j] = std::get<0>(_e0[vertex_offsets[i] + j]);
      }
    });
    vertex_offsets[n] = v_data[n - 1].offset + v_data[n - 1].degree;
    if constexpr (shuffle) {
      gbbs::parallel_for(0, n, [&](size_t i) {
        std::mt19937 g(i);
        std::shuffle(edges.data() + vertex_offsets[i],
                     edges.data() + vertex_offsets[i + 1], g);
      });
    }
  }

  // Move constructor
  CSR(CSR &&other) noexcept {
    n = other.n;
    m = other.m;
    other.n = 0;
    other.m = 0;
    vertex_offsets = other.vertex_offsets;
    vertex_weights = other.vertex_weights;
    other.vertex_weights = nullptr;
    deletion_fn = std::move(other.deletion_fn);
    edges = other.edges;
  }

  // Move assignment
  CSR &operator=(CSR &&other) noexcept {
    n = other.n;
    m = other.m;
    other.n = 0;
    other.m = 0;
    vertex_offsets = other.vertex_offsets;
    vertex_weights = other.vertex_weights;
    other.vertex_weights = nullptr;
    deletion_fn = std::move(other.deletion_fn);
    edges = other.edges;
  }

  // Copy constructor
  CSR(const CSR &other) : n(other.n), m(other.m) {
    debug(std::cout << "Copying symmetric graph." << std::endl;);
    vertex_offsets.reset((edge_t *)malloc((n + 1) * sizeof(edge_t)));
    edges.reset((node_t *)malloc((m) * sizeof(node_t)));
    vertex_weights = nullptr;
    if (other.vertex_weights != nullptr) {
      vertex_weights = gbbs::new_array_no_init<vertex_weight_type>(n);
      parallel_for(
          0, n, [&](size_t i) { vertex_weights[i] = other.vertex_weights[i]; });
    }
    deletion_fn = [&]() {
      if (vertex_weights != nullptr) {
        gbbs::free_array(vertex_weights, n);
      }
    };
    parallel_for(0, n + 1, [&](size_t i) {
      vertex_offsets[i] = other.vertex_offsets[i];
    });
    parallel_for(0, m, [&](size_t i) { edges[i] = other.edges[i]; });
  }

  size_t get_memory_size() {
    return (n + 1) * sizeof(edge_t) + m * sizeof(node_t) + sizeof(*this);
  }
  ~CSR() { deletion_fn(); }
  size_t N() const { return n; }
  size_t M() const { return m; }

private:
  struct free_delete {
    void operator()(void *x) { free(x); }
  };
  std::unique_ptr<edge_t[], free_delete> vertex_offsets = nullptr;
  std::unique_ptr<node_t[], free_delete> edges = nullptr;
  node_t n = 0;
  edge_t m = 0;
  // called to delete the graph
  vertex_weight_type *vertex_weights = nullptr;
  std::function<void()> deletion_fn = []() {};
  // Pointer to vertex weights
};

#ifdef LONG
using edge_t = uint64_t;
#else
using edge_t = uint32_t;
#endif
#ifdef EDGELONG
using node_t = uint64_t;
#else
using node_t = uint32_t;
#endif
using graph_impl = CSR<node_t, edge_t>;

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
    std::cerr << "reads in uncompressed files\n";
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