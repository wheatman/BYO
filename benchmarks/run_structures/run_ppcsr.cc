

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

#include <cstddef>
#include <functional>

#ifdef HOMEGROWN
#define PARLAY 1
#endif

#include "gbbs/bridge.h"
#include "terrace/PMA.hpp"

#include "../run_unweighted.h"

template <class W> struct symmetric_ppcsr_graph {
  using weight_type = W;
  static constexpr bool binary = true;
  static_assert(binary);
  using vertex_weight_type = double;
  using edge_type = std::tuple<gbbs::uintE, W>;

  // size_t M() const { return nodes.get_num_edges(); }

  size_t N() const { return nodes.nodes.size(); }

  auto degree(size_t i) const { return nodes.nodes[i].num_neighbors; }

  template <class F> void map_neighbors(size_t i, F f) const {
    uint64_t start = nodes.nodes[i].beginning + 1;
    uint64_t end = nodes.nodes[i].end;
    W empty_weight = W();
    for (uint64_t j = start; j < end; j++) {
      if (nodes.edges.dests[j] != NULL_VAL) {
        f(i, nodes.edges.dests[j], empty_weight);
      }
    }
  }
  template <class F> void map_neighbors_early_exit(size_t i, F f) const {
    uint64_t start = nodes.nodes[i].beginning + 1;
    uint64_t end = nodes.nodes[i].end;
    W empty_weight = W();
    for (uint64_t j = start; j < end; j++) {
      if (nodes.edges.dests[j] != NULL_VAL) {
        if (f(i, nodes.edges.dests[j], empty_weight)) {
          break;
        }
      }
    }
  }
  template <class F> void parallel_map_neighbors(size_t i, F f) const {
    uint64_t start = nodes.nodes[i].beginning + 1;
    uint64_t end = nodes.nodes[i].end;
    W empty_weight = W();
    gbbs::parallel_for(start, end, [&](auto j) {
      if (nodes.edges.dests[j] != NULL_VAL) {
        f(i, nodes.edges.dests[j], empty_weight);
      }
    });
  }

  template <class F, typename Block_F = std::nullptr_t>
  void parallel_map_neighbors_early_exit(size_t i, F f,
                                         Block_F block_check = {}) const {
    uint64_t start = nodes.nodes[i].beginning + 1;
    uint64_t end = nodes.nodes[i].end;
    uint64_t size = end - start;
    if (size < 1000) {
      map_neighbors(i, f);
    }
    size_t b_size = 2048;
    size_t n_blocks = size / b_size + 1;
    weight_type empty_weight{};
    gbbs::parallel_for(0, n_blocks, [&](size_t b) {
      if constexpr (!std::is_same_v<Block_F, std::nullptr_t>) {
        if (!block_check()) {
          return;
        }
      }
      size_t range_start = start + b * b_size;
      size_t range_end =
          start + std::min((b + 1) * b_size, static_cast<size_t>(size));
      for (size_t j = range_start; j < range_end; j++) {
        if (nodes.edges.dests[j] != NULL_VAL) {
          if (f(i, nodes.edges.dests[j], empty_weight)) {
            break;
          }
        }
      }
    });
  }
  // ======================= Constructors and fields  ========================
  symmetric_ppcsr_graph() : vertex_weights(nullptr) {}
  // for now use the same constructor as the exsisting one for ease, but
  // probably write a more general constuctor for everyone to use later
  // TODO(wheatman) figure out what to do with _deletion_fn, probably should
  // just use unique pointer
  symmetric_ppcsr_graph(auto *v_data, size_t n, size_t m,
                        std::function<void()> _deletion_fn, edge_type *_e0,
                        vertex_weight_type *_vertex_weights = nullptr)
      : nodes(n), vertex_weights(_vertex_weights), deletion_fn(_deletion_fn) {

    gbbs::sequence<uint32_t> srcs = gbbs::sequence<uint32_t>::uninitialized(m);
    gbbs::sequence<uint32_t> dests = gbbs::sequence<uint32_t>::uninitialized(m);
    gbbs::sequence<uint32_t> degrees(n, 0);
    gbbs::sequence<uint8_t> all_true(m, 1);

    printf("convert to edge list\n");
    gbbs::parallel_for(0, n, [&](uint32_t i) {
      for (size_t j = v_data[i].offset; j < v_data[i].offset + v_data[i].degree;
           j++) {
        srcs[j] = i;
        dests[j] = std::get<0>(_e0[j]);
        // TODO: weights if necessary
      }
    });
    nodes.build_from_edges(srcs.data(), dests.data(), all_true.data(), n, m,
                           degrees.data());

    // std::cout << "verifying edges\n";

    // for (size_t i = 0; i < n; i++) {
    //   size_t edge_offset = v_data[i].offset;
    //   if (degree(i) != v_data[i].degree) {
    //     std::cout << "bad degree for node " << i << " got " << degree(i)
    //               << " expected " << v_data[i].degree << "\n";
    //   }

    //   map_neighbors(i, [&](auto src, auto dest, auto val) {
    //     if (src != i) {
    //       std::cout << "1. bad edge (" << src << ", " << dest << ")\n";
    //     }
    //     if (dest != std::get<0>(_e0[edge_offset])) {
    //       std::cout << "2. bad edge (" << src << ", " << dest << ")\n";
    //       std::cout << "the " << edge_offset - v_data[i].offset
    //                 << " edge of the node\n";
    //       std::cout << "the node has degree " << degree(i) << "\n";
    //       std::cout << "edge_offset = " << edge_offset << "\n";
    //     }
    //     edge_offset += 1;
    //   });
    // }

    // printf("starting new build from batch\n");
    // nodes.build_from_batch(srcs.data(), dests.data(), n, m);
    // // nodes.add_edge_batch_no_perm(srcs.data(), dests.data(), m);

    // printf("** STARTING VERIFY **\n");
    // for (uint32_t i = 0; i < n; i++) {
    //   nodes.verify_neighbors(i, dests.data() + v_data[i].offset);
    // }
    // printf("** FINISHED VERIFY **\n");

    // printf("AFTER INIT nodes in graph = %u\n", nodes.get_num_vertices());
    // printf("AFTER INIT edges in graph = %lu\n", nodes.get_num_edges());
  }

  // TODO: write these
  // Move constructor
  // symmetric_terrace_graph(symmetric_terrace_graph &&other) noexcept {
  //   printf("move constructor\n");
  //   /*
  //     nodes = other.nodes;
  //     vertex_weights = other.vertex_weights;
  //     other.vertex_weights = nullptr;
  //   */
  // }
  // /*
  // // Move assignment
  // symmetric_terrace_graph &
  // operator=(symmetric_terrace_graph &&other) noexcept {
  //   nodes = other.nodes;
  //   vertex_weights = other.vertex_weights;
  //   other.vertex_weights = nullptr;
  // }
  // */
  // // Copy constructor
  // symmetric_terrace_graph(const symmetric_terrace_graph &other)
  //     : nodes(other.N()) {
  //   printf("copy constructor\n");
  //   for (size_t i = 0; i < N(); i++) {
  //     other.map_neighbors(
  //         i, [&](auto src, auto dest, auto val) { nodes.add_edge(src, dest);
  //         });
  //   }

  //   vertex_weights = nullptr;
  //   if (other.vertex_weights != nullptr) {
  //     vertex_weights = gbbs::new_array_no_init<vertex_weight_type>(N());
  //     gbbs::parallel_for(0, N(), [&](size_t i) {
  //       vertex_weights[i] = other.vertex_weights[i];
  //     });
  //   }
  //   deletion_fn = [&]() {
  //     if (vertex_weights != nullptr) {
  //       gbbs::free_array(vertex_weights, N());
  //     }
  //   };
  // }
  ~symmetric_ppcsr_graph() { deletion_fn(); }

  void insert_batch(std::tuple<uint32_t, uint32_t> *es, size_t n) {
    nodes.add_edge_batch_wrapper(std::bit_cast<pair_uint *>(es), n);
  }

  void remove_batch(std::tuple<uint32_t, uint32_t> *es, size_t n) {
    nodes.remove_edge_batch_wrapper(std::bit_cast<pair_uint *>(es), n);
  }

  // Graph Data
  graphstore::PMA nodes;
  vertex_weight_type *vertex_weights;
  // called to delete the graph
  std::function<void()> deletion_fn;
};

using graph_impl = symmetric_ppcsr_graph<gbbs::empty>;

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
