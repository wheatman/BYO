

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

#include <cstdint>
#include <functional>

#ifdef HOMEGROWN
#define PARLAY 1
#endif

#include "reducer.hpp"
#define NO_TLX
#include "PMA/PCSR.hpp"
#include "gbbs/bridge.h"

template <template <class W> class vertex_type, class W>
struct symmetric_PCSR_graph {
  using vertex = vertex_type<W>;
  using weight_type = W;
  static constexpr bool binary = std::is_same_v<gbbs::empty, W>;
  static_assert(binary);
  using vertex_weight_type = double;
  using edge_type = typename vertex::edge_type;

  size_t N() const { return nodes.num_nodes(); }

  auto degree(size_t i) const { return nodes.get_degree(i); }

  template <class F> void map_neighbors(size_t i, F f) const {
    W empty_weight = W();
    nodes.map_neighbors<-1>(
        i, [&](auto src, auto dest) { return f(src, dest, empty_weight); },
        nullptr, false);
  }

  template <class F> void parallel_map_neighbors(size_t i, F f) const {
    W empty_weight = W();
    nodes.map_neighbors<-1>(
        i, [&](auto src, auto dest) { return f(src, dest, empty_weight); },
        nullptr, true);
  }

  template <class F> void map_neighbors_early_exit(size_t i, F f) const {
    W empty_weight = W();
    nodes.map_neighbors<1>(
        i, [&](auto src, auto dest) { return f(src, dest, empty_weight); },
        nullptr, false);
  }

  template <class F>
  void parallel_map_neighbors_early_exit(size_t i, F f) const {
    W empty_weight = W();
    nodes.map_neighbors<1>(
        i, [&](auto src, auto dest) { return f(src, dest, empty_weight); },
        nullptr, true);
  }

  // ======================= Constructors and fields  ========================
  symmetric_PCSR_graph() : vertex_weights(nullptr) {}

  // for now use the same constructor as the exsisting one for ease, but
  // probably write a more general constuctor for everyone to use later
  // TODO(wheatman) figure out what to do with _deletion_fn, probably should
  // just use unique pointer
  symmetric_PCSR_graph(auto *v_data, size_t n, size_t m,
                       std::function<void()> _deletion_fn, edge_type *_e0,
                       vertex_weight_type *_vertex_weights = nullptr)
      : nodes(n), vertex_weights(_vertex_weights), deletion_fn(_deletion_fn) {

    Reducer_Vector<std::tuple<uint32_t, uint32_t>> vec;
    gbbs::parallel_for(0, n, [&](uint64_t i) {
      for (size_t j = v_data[i].offset; j < v_data[i].offset + v_data[i].degree;
           j++) {
        vec.push_back({i, std::get<0>(_e0[j])});
      }
    });
    auto seq = vec.get_sequence();
    nodes.insert_batch(seq);
  }

  void insert_sorted_batch(auto *es, size_t n) {
    auto slice = parlay::slice(es, es + n);
    nodes.insert_batch(slice, true);
  }

  void remove_sorted_batch(auto *es, size_t n) {
    auto slice = parlay::slice(es, es + n);
    nodes.remove_batch(slice, true);
  }

  // Move constructor
  symmetric_PCSR_graph(symmetric_PCSR_graph &&other) noexcept {
    nodes = other.nodes;
    vertex_weights = other.vertex_weights;
    other.vertex_weights = nullptr;
    deletion_fn = std::move(other.deletion_fn);
  }

  // Move assignment
  symmetric_PCSR_graph &operator=(symmetric_PCSR_graph &&other) noexcept {
    nodes = other.nodes;
    vertex_weights = other.vertex_weights;
    other.vertex_weights = nullptr;
    deletion_fn = std::move(other.deletion_fn);
  }

  // Copy constructor
  symmetric_PCSR_graph(const symmetric_PCSR_graph &other) : nodes(other.nodes) {
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

  size_t get_memory_size() { return nodes.get_memory_size(); }

  ~symmetric_PCSR_graph() { deletion_fn(); }

  // Graph Data
  using traits = PMA_traits<uncompressed_leaf<uint32_t>, Eytzinger, 0, false,
                            false, false, 0, true, true>;
  PCSR<traits> nodes;
  vertex_weight_type *vertex_weights;
  // called to delete the graph
  std::function<void()> deletion_fn;
};

#include "../run_unweighted.h"

using graph_impl = symmetric_PCSR_graph<gbbs::symmetric_vertex, gbbs::empty>;

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
    std::cerr << "is not compressed, and reads in uncompressed files\n";
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