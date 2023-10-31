

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

#include <cstdio>
#include <cstdint>
#include <iterator>
#include <utility>

#define PARLAYSCHE 1
#include "gbbs/macros.h"
#include "graph/api.h"
#include "pbbslib/get_time.h"

class AspenWrapper {
  /*
  // using list = compressed_lists;
  using AT = compressed_lists::array_type;

  struct edge_entry {
    using key_t = uintV; // the 'head' edge of this node.
    using val_t = AT*; // the set of edges stored in this node.
    static bool comp(const key_t& a, const key_t& b) { return a < b; }
    using aug_t = uintV; // num. edges in this subtree
    static aug_t get_empty() { return 0; }
    static aug_t from_entry(const key_t& k, const val_t& v) { return 1 + compressed_lists::node_size(v); }
    static aug_t combine(const aug_t& a, const aug_t& b) { return a + b; }
    using entry_t = std::pair<key_t,val_t>;
    static entry_t copy_entry(const entry_t& e) {
      // TODO: Instead of copying, bump a ref-ct (note that copy_node and
      // deallocate can implement these semantics internally)
      return std::make_pair(e.first, compressed_lists::copy_node(e.second));
    }
    static void del(entry_t& e) {
      if (e.second) {
        // TODO: Should decrement ref-ct, free if only owner
        compressed_lists::deallocate(e.second);
      }
    }
  };
  */

  // using edge_list = aug_map<edge_entry>;
  using edge_tree = tree_plus::treeplus;
  edge_tree edges;

  struct aux_init{
    struct singleton{
      singleton(){
        using edge_list = tree_plus::edge_list;
        edge_list::init();
        // following the setting in the original Aspen code
        const size_t n_max = 70000000;
        edge_list::reserve(n_max/16);
        lists::init(n_max);
      }
    };
    
    aux_init() { static singleton impl; }
  };

  // make sure edge_list always gets initialized if any edge_tree is used
  aux_init aux = {};

 public:
  AspenWrapper() = default;
  AspenWrapper(const AspenWrapper&) = default;
  AspenWrapper(AspenWrapper&&) = default;
  AspenWrapper& operator=(const AspenWrapper&) = default;
  AspenWrapper& operator=(AspenWrapper&&) = default;

  template<typename Iter>
  AspenWrapper(Iter begin, Iter end) {
    edges.del();
    edges = edge_tree(parlay::slice(begin,end), 0);
  }

  size_t size() const { return edges.size(); }

  template <class F>
  void map(F f) const {
    auto g = [&](const uintV& ngh) { f(ngh, {}); };
    edges.iter_elms(0, g);
  }

  template <class F>
  void map_early_exit(F f) const {
    auto g = [&](const uintV& ngh) { return f(ngh, {}); };
    edges.iter_elms_cond(0, g);
  }

  template <class F>
  void parallel_map(F f) const {
    auto g = [&](const uintV& ngh, size_t ind) { f(ngh, {}); };
    edges.map_elms(0, g);
  }

  //template <class F>
  //void parallel_map_early_exit(F f) const {
    //auto map_f = [&](const auto& et) { return f(et, {}); };
    //edges.foreach_cond_par(edges, map_f, []() { return true; });
  //}

  void insert_sorted_batch(auto es, size_t n) {
    edges = tree_plus::uniont(edges, edge_tree(parlay::slice(es,es+n),0), 0);
  }
  void remove_sorted_batch(auto es, size_t n) {
    // auto to_remove = parlay::tabulate(n, [&](size_t i){return *(es+i);});
    // auto bool_seq = parlay::delayed_seq<bool>(to_remove.size(), [&](size_t i) {
    //   return (i == 0 || to_remove[i] != to_remove[i-1]);
    // });
    // to_remove = parlay::pack(to_remove, bool_seq);
    edges = tree_plus::difference(edge_tree(parlay::slice(es,es+n),0), edges, 0);
  }
};

#include "../run_unweighted.h"

#ifdef USE_INPLACE
static constexpr bool use_inplace = true;
#else
static constexpr bool use_inplace = false;
#endif

using graph_impl = gbbs::graph_implementations::symmetric_set_graph<
    gbbs::symmetric_vertex, gbbs::empty, AspenWrapper,
    /* inplace */ use_inplace>;

using graph_api = gbbs::full_api;

using graph_t = gbbs::Graph<graph_impl, /* symmetric */ true, graph_api>;

int main(int argc, char* argv[]) {
  gbbs::commandLine P(argc, argv, " [-s] <inFile>");
  char* iFile = P.getArgument(0);
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
    std::cerr << "is always compressed, but reads in uncompressed files\n";
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
