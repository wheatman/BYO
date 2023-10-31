

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

#define NO_TLX
#include "PMA/CPMA.hpp"
#include "gbbs/bridge.h"
#include "reducer.hpp"

namespace {
template <class F, bool no_early_exit_, class W> struct F2 {
  F f;
  F2(F f_) : f(f_) {}
  W empty_weight = W();
  auto operator()(auto src, auto dest) { return f(src, dest, empty_weight); }
  static constexpr bool no_early_exit = no_early_exit_;
};
}; // namespace
template <template <class W> class vertex_type, class W> struct PMA_graph {
  using vertex = vertex_type<W>;
  using vertex_weight_type = double;
  using weight_type = W;
  using edge_type = typename vertex::edge_type;

  size_t N() const { return pma.num_nodes(); }

  template <class F> void map_neighbors(size_t i, F f) const {
    pma.map_neighbors(i, F2<F, true, W>(f), extra_data, false);
  }

  template <class F> void parallel_map_neighbors(size_t i, F f) const {
    pma.map_neighbors(i, F2<F, true, W>(f), extra_data, true);
  }

  template <class F> void map_neighbors_early_exit(size_t i, F f) const {
    pma.map_neighbors(i, F2<F, false, W>(f), extra_data, false);
  }

  template <class F>
  void parallel_map_neighbors_early_exit(size_t i, F f) const {
    pma.map_neighbors(i, F2<F, false, W>(f), extra_data, true);
  }

  void insert_sorted_batch(std::tuple<uint32_t, uint32_t> *es, size_t n) {
    pma.insert_batch(MultiPointer<uint64_t>(std::bit_cast<uint64_t *>(es)), n, true);
  }

  void remove_sorted_batch(std::tuple<uint32_t, uint32_t> *es, size_t n) {
    pma.remove_batch(std::bit_cast<uint64_t *>(es), n, true);
  }

  // ======================= Constructors and fields  ========================
  PMA_graph() : vertex_weights(nullptr) {}

  PMA_graph(auto *v_data, size_t n, size_t m,
            std::function<void()> _deletion_fn, edge_type *_e0,
            vertex_weight_type *_vertex_weights = nullptr)
      : vertex_weights(_vertex_weights), deletion_fn(_deletion_fn) {
    Reducer_Vector<uint64_t> vec;
    gbbs::parallel_for(0, n, [&](uint64_t i) {
      for (size_t j = v_data[i].offset; j < v_data[i].offset + v_data[i].degree;
           j++) {
        vec.push_back(i << 32UL | std::get<0>(_e0[j]));
      }
    });
    auto seq = vec.get_sequence();
    pma.insert_batch(seq.data(), seq.size());
    extra_data = pma.getExtraData();
  }

  // Move constructor
  PMA_graph(PMA_graph &&other) noexcept {
    pma = other.pma;
    vertex_weights = other.vertex_weights;
    other.vertex_weights = nullptr;
    extra_data = other.extra_data;
    deletion_fn = std::move(other.deletion_fn);
  }

  // Move assignment
  PMA_graph &operator=(PMA_graph &&other) noexcept {
    pma = other.pma;
    vertex_weights = other.vertex_weights;
    other.vertex_weights = nullptr;
    extra_data = other.extra_data;
    deletion_fn = std::move(other.deletion_fn);
  }

  // Copy constructor
  PMA_graph(const PMA_graph &other) : pma(other.pma) {
    vertex_weights = nullptr;
    if (other.vertex_weights != nullptr) {
      vertex_weights = gbbs::new_array_no_init<vertex_weight_type>(N());
      parallel_for(0, N(), [&](size_t i) {
        vertex_weights[i] = other.vertex_weights[i];
      });
    }
    extra_data = pma.getExtraData();
    deletion_fn = [&]() {
      if (vertex_weights != nullptr) {
        gbbs::free_array(vertex_weights, N());
      }
    };
  }

  ~PMA_graph() { deletion_fn(); }

  // Graph Data
  CPMA<PMA_SETTINGS<uint64_t>> pma;
  vertex_weight_type *vertex_weights;
  decltype(pma.getExtraData()) extra_data;
  // called to delete the graph
  std::function<void()> deletion_fn;
};

#include "../run_unweighted.h"

using graph_impl = PMA_graph<gbbs::symmetric_vertex, gbbs::empty>;

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
    std::cerr << "reads in uncompressed files\n";
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