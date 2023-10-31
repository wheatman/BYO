

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

#include "dhb/dynamic_hashed_blocks.h"
#include "gbbs/bridge.h"
#include "semisort.h"

template <template <class W> class vertex_type, class W>
struct symmetric_dhb_graph {
  using vertex = vertex_type<W>;
  using weight_type = W;
  static constexpr bool binary = std::is_same_v<gbbs::empty, W>;
  using vertex_weight_type = double;
  using edge_type = typename vertex::edge_type;

  size_t N() const { return nodes.vertices_count(); }

  auto degree(size_t i) const { return nodes.degree(i); }

  template <class F> void map_neighbors(size_t i, F f) const {
    auto neighbors = nodes.neighbors(i);
    for (auto neighbor : neighbors) {
      size_t other_vertex = neighbor.vertex();
      weight_type weight = neighbor.data();
      f(i, other_vertex, weight);
    }
  }

  template <class F> void map_neighbors_early_exit(size_t i, F f) const {
    auto neighbors = nodes.neighbors(i);
    for (auto neighbor : neighbors) {
      size_t other_vertex = neighbor.vertex();
      weight_type weight = neighbor.data();
      if (f(i, other_vertex, weight)) {
        return;
      }
    }
  }
  // ======================= Constructors and fields  ========================
  symmetric_dhb_graph() : vertex_weights(nullptr) {}

  // for now use the same constructor as the exsisting one for ease, but
  // probably write a more general constuctor for everyone to use later
  // TODO(wheatman) figure out what to do with _deletion_fn, probably should
  // just use unique pointer
  symmetric_dhb_graph(auto *v_data, size_t n, size_t m,
                      std::function<void()> _deletion_fn, edge_type *_e0,
                      vertex_weight_type *_vertex_weights = nullptr)
      : nodes(n), vertex_weights(_vertex_weights), deletion_fn(_deletion_fn) {

    gbbs::parallel_for(0, n, [&](uint32_t i) {
      for (size_t j = v_data[i].offset; j < v_data[i].offset + v_data[i].degree;
           j++) {
        if constexpr (binary) {
          nodes.insert(i, std::get<0>(_e0[j]), {}, true);
        } else {
          nodes.insert(i, std::get<0>(_e0[j]), std::get<1>(_e0[j]), true);
        }
      }
    });
  }

  // Move constructor
  symmetric_dhb_graph(symmetric_dhb_graph &&other) noexcept {
    nodes = other.nodes;
    vertex_weights = other.vertex_weights;
    other.vertex_weights = nullptr;
    deletion_fn = std::move(other.deletion_fn);
  }

  // Move assignment
  symmetric_dhb_graph &operator=(symmetric_dhb_graph &&other) noexcept {
    nodes = other.nodes;
    vertex_weights = other.vertex_weights;
    other.vertex_weights = nullptr;
    deletion_fn = std::move(other.deletion_fn);
  }

  // Copy constructor
  symmetric_dhb_graph(const symmetric_dhb_graph &other) : nodes(other.nodes) {
    vertex_weights = nullptr;
    if (other.vertex_weights != nullptr) {
      vertex_weights = gbbs::new_array_no_init<vertex_weight_type>(N());
      parallel_for(0, N(), [&](size_t i) {
        vertex_weights[i] = other.vertex_weights[i];
      });
    }
    deletion_fn = [&]() {
      if (vertex_weights != nullptr) {
        gbbs::free_array(vertex_weights, nodes.size());
      }
    };
  }

  void insert_sorted_grouped_batch(const parlay::sequence<std::pair<uint32_t, parlay::sequence<uint32_t>>> & els) {
    parlay::parallel_for(0, els.size(), [&](size_t i) {
      for (const auto& el : els[i].second) {
          nodes.insert(els[i].first, el, {}, true);
      }
    });
  }

  void remove_sorted_grouped_batch(const parlay::sequence<std::pair<uint32_t, parlay::sequence<uint32_t>>> & els) {
    parlay::parallel_for(0, els.size(), [&](size_t i) {
      for (const auto& el : els[i].second) {
          nodes.removeEdge(els[i].first, el);
      }
    });
  }

  ~symmetric_dhb_graph() { deletion_fn(); }

  // Graph Data
  dhb::Matrix<weight_type> nodes;
  vertex_weight_type *vertex_weights;
  // called to delete the graph
  std::function<void()> deletion_fn;
};

#include "../run_unweighted.h"

using graph_impl = symmetric_dhb_graph<gbbs::symmetric_vertex, gbbs::empty>;

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