// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <cstddef>
#include <stdlib.h>
#include <fstream>
#include <functional>
#include <iostream>
#include <optional>
#include <string>
#include <type_traits>

#include <vector>
#include <functional>
#include <random>
#include <algorithm>
#include <ranges>

#include "bridge.h"
#include "compressed_vertex.h"
#include "edge_array.h"
#include "flags.h"
#include "macros.h"
#include "vertex.h"

#include "reducer.hpp"

#include "semisort.h"

namespace gbbs {

  namespace graph_implementations {

//  Compressed Sparse Row (CSR) based representation for symmetric graphs.
//  Takes two template parameters:
//  1) vertex_type: vertex template, parametrized by the weight type associated
//  with each edge
//  2) W: the edge weight template
//  The graph is represented as an array of edges of type
//  vertex_type::edge_type.
//  For uncompressed vertices, this type is equal to tuple<uintE, W>.
template <template <class W> class vertex_type, class W, bool shuffle = false>
struct symmetric_graph {
  using vertex = vertex_type<W>;
  using weight_type = W;
  using edge_type = typename vertex::edge_type;
  using graph = symmetric_graph<vertex_type, W>;
  using vertex_weight_type = double;

  size_t num_vertices() const { return n; }
  size_t num_edges() const { return m; }

  // ======== Graph operators that perform packing ========
  template <class P>
  uintE packNeighbors(uintE id, P& p, uint8_t* tmp) {
    uintE new_degree =
        get_vertex(id).out_neighbors().pack(p, (std::tuple<uintE, W>*)tmp);
    v_data[id].degree = new_degree;  // updates the degree
    return new_degree;
  }

  // degree must be <= old_degree
  void decreaseVertexDegree(uintE id, uintE degree) {
    assert(degree <= v_data[id].degree);
    v_data[id].degree = degree;
  }

  void zeroVertexDegree(uintE id) { decreaseVertexDegree(id, 0); }

  sequence<std::tuple<uintE, uintE, W>> edges() const {
    using g_edge = std::tuple<uintE, uintE, W>;
    auto degs = sequence<size_t>::from_function(
        n, [&](size_t i) { return get_vertex(i).out_degree(); });
    size_t sum_degs = parlay::scan_inplace(make_slice(degs));
    assert(sum_degs == m);
    auto edges = sequence<g_edge>(sum_degs);
    parallel_for(0, n,
                 [&](size_t i) {
                   size_t k = degs[i];
                   auto map_f = [&](const uintE& u, const uintE& v,
                                    const W& wgh) {
                     edges[k++] = std::make_tuple(u, v, wgh);
                   };
                   get_vertex(i).out_neighbors().map(map_f, false);
                 },
                 1);
    return edges;
  }

  template <class F>
  void mapEdges(F f, bool parallel_inner_map = true, size_t granularity = 1) const {
    parallel_for(0, n,
                 [&](size_t i) {
                   get_vertex(i).out_neighbors().map(f, parallel_inner_map);
                 },
                 granularity);
  }

  template <class F>
  void map_neighbors(size_t i, F f) const {
    get_vertex(i).out_neighbors().map(f, false);
  }
  template <class F>
  void map_neighbors_early_exit(size_t i, F f) const {
    constexpr bool support_map_early_exit =
        requires() {
      get_vertex(i).out_neighbors().map_early_exit(f, false);
    };
    if constexpr(support_map_early_exit) {
      get_vertex(i).out_neighbors().map_early_exit(f, false);
    } else {
      auto neighbors = get_vertex(i).out_neighbors();
      for (size_t j = 0; j < neighbors.degree; j++) {
        auto target = neighbors.get_neighbor(j);
        auto weight = neighbors.get_weight(j);
        if (f(i, target, weight)){
          break;
        }
      }
    }
  }
  template <class F>
  void parallel_map_neighbors(size_t i, F f) const {
    get_vertex(i).out_neighbors().map(f, true);
  }

  template <class F, typename Block_F = std::nullptr_t>
  void parallel_map_neighbors_early_exit(size_t i, F f, Block_F block_check = {}) const {
        constexpr bool support_map_early_exit =
        requires() {
      get_vertex(i).out_neighbors().map_early_exit(f, true);
    };
    if constexpr(support_map_early_exit) {
      get_vertex(i).out_neighbors().map_early_exit(f, true);
    } else {
      auto neighbors = get_vertex(i).out_neighbors();
      if (neighbors.degree < 1000) {
        map_neighbors_early_exit(i, f);
        return;
      }
      size_t b_size = 2048;
      size_t n_blocks = neighbors.degree / b_size + 1;
      parallel_for (0, n_blocks, [&] (size_t b){
        if constexpr(!std::is_same_v<Block_F, std::nullptr_t>) {
          if (!block_check()) {
            return;
          }
        }
        size_t start = b * b_size;
        size_t end = std::min((b + 1) * b_size,
                        static_cast<size_t>(neighbors.degree));
        for (size_t j = start; j < end; j++) {
          auto target = neighbors.get_neighbor(j);
          auto weight = neighbors.get_weight(j);
          if (f(i, target, weight)){
            return;
          }
        }
      });
    }
  }

  uintE degree(size_t i) const {
    return get_vertex(i).out_neighbors().degree;
  }

  template <class M, class R>
  typename R::T reduceEdges(M map_f, R reduce_f) const {
    using T = typename R::T;
    auto D = parlay::delayed_seq<T>(n, [&](size_t i) {
      return get_vertex(i).out_neighbors().reduce(map_f, reduce_f);
    });
    return parlay::reduce(D, reduce_f);
  }

  // ======================= Constructors and fields  ========================
  symmetric_graph()
      : v_data(nullptr),
        e0(nullptr),
        vertex_weights(nullptr),
        n(0),
        m(0),
        deletion_fn([]() {}) {}

  symmetric_graph(vertex_data *v_data, size_t n, size_t m,
                  std::function<void()> _deletion_fn, edge_type *_e0,
                  vertex_weight_type *_vertex_weights = nullptr)
      : v_data(v_data), e0(_e0), vertex_weights(_vertex_weights), n(n), m(m),
        deletion_fn(_deletion_fn) {
    if constexpr (shuffle) {
      auto degs = sequence<size_t>::from_function(
          n, [&](size_t i) { return get_vertex(i).out_degree(); });
      size_t sum_degs = parlay::scan_inplace(make_slice(degs));
      assert(sum_degs == m);
      parallel_for(
          0, n,
          [&](size_t i) {
            size_t start = degs[i];
            size_t end = sum_degs;
            if (i < n - 1) {
              end = degs[i + 1];
            }
            std::mt19937 g(i);
            std::shuffle(e0 + start, e0 + end, g);
          },
          1);
    }
  }

  // Move constructor
  symmetric_graph(symmetric_graph&& other) noexcept {
    n = other.n;
    m = other.m;
    v_data = other.v_data;
    e0 = other.e0;
    vertex_weights = other.vertex_weights;
    deletion_fn = std::move(other.deletion_fn);
    other.v_data = nullptr;
    other.e0 = nullptr;
    other.vertex_weights = nullptr;
    other.deletion_fn = []() {};
  }

  // Move assignment
  symmetric_graph& operator=(symmetric_graph&& other) noexcept {
    n = other.n;
    m = other.m;
    v_data = other.v_data;
    e0 = other.e0;
    vertex_weights = other.vertex_weights;
    deletion_fn = std::move(other.deletion_fn);
    other.v_data = nullptr;
    other.e0 = nullptr;
    other.vertex_weights = nullptr;
    other.deletion_fn = []() {};
    return *this;
  }

  // Copy constructor
  symmetric_graph(const symmetric_graph& other) {
    n = other.n;
    m = other.m;
    v_data = gbbs::new_array_no_init<vertex_data>(n);
    e0 = gbbs::new_array_no_init<edge_type>(m);
    parallel_for(0, n, [&](size_t i) { v_data[i] = other.v_data[i]; });
    parallel_for(0, m, [&](size_t i) { e0[i] = other.e0[i]; });
    deletion_fn = [&]() {
      gbbs::free_array(v_data, n);
      gbbs::free_array(e0, m);
      if (vertex_weights != nullptr) {
        gbbs::free_array(vertex_weights, n);
      }
    };
    vertex_weights = nullptr;
    if (other.vertex_weights != nullptr) {
      vertex_weights = gbbs::new_array_no_init<vertex_weight_type>(n);
      parallel_for(
          0, n, [&](size_t i) { vertex_weights[i] = other.vertex_weights[i]; });
    }
  }

  ~symmetric_graph() { deletion_fn(); }

  vertex get_vertex(uintE i) const { return vertex(e0, v_data[i], i); }

  // Graph Data
  vertex_data* v_data;
  // Pointer to edges
  edge_type* e0;
  // Pointer to vertex weights
  vertex_weight_type* vertex_weights;

  // number of vertices in G
  size_t n;
  // number of edges in G
  size_t m;

  // called to delete the graph
  std::function<void()> deletion_fn;
};

// Similar to symmetric_graph, but edges are not necessarily allocated
// consecutively. The structure simply stores an array of vertex
// objects (which store an 8-byte pointer, and a uintE degree each).
template <template <class W> class vertex_type, class W>
struct symmetric_ptr_graph {
  using vertex = vertex_type<W>;
  using weight_type = W;
  using edge_type = typename vertex::edge_type;
  using graph = symmetric_ptr_graph<vertex_type, W>;
  using vertex_weight_type = double;

  size_t num_vertices() const { return n; }
  size_t num_edges() const { return m; }

  // ======== Graph operators that perform packing ========
  template <class P>
  uintE packNeighbors(uintE id, P& p, uint8_t* tmp) {
    uintE new_degree =
        get_vertex(id).out_neighbors().pack(p, (std::tuple<uintE, W>*)tmp);
    vertices[id].degree = new_degree;  // updates the degree
    return new_degree;
  }

  // degree must be <= old_degree
  void decreaseVertexDegree(uintE id, uintE degree) {
    assert(degree <= vertices[id].degree);
    vertices[id].degree = degree;
  }

  void zeroVertexDegree(uintE id) { decreaseVertexDegree(id, 0); }

  sequence<std::tuple<uintE, uintE, W>> edges() const {
    using g_edge = std::tuple<uintE, uintE, W>;
    auto degs = sequence<size_t>::from_function(
        n, [&](size_t i) { return get_vertex(i).out_degree(); });
    size_t sum_degs = parlay::scan_inplace(make_slice(degs));
    assert(sum_degs == m);
    auto edges = sequence<g_edge>(sum_degs);
    parallel_for(0, n,
                 [&](size_t i) {
                   size_t k = degs[i];
                   auto map_f = [&](const uintE& u, const uintE& v,
                                    const W& wgh) {
                     edges[k++] = std::make_tuple(u, v, wgh);
                   };
                   get_vertex(i).out_neighbors().map(map_f, false);
                 },
                 1);
    return edges;
  }

  template <class F>
  void mapEdges(F f, bool parallel_inner_map = true) const {
    parallel_for(0, n,
                 [&](size_t i) {
                   get_vertex(i).out_neighbors().map(f, parallel_inner_map);
                 },
                 1);
  }

  template <class M, class R>
  typename R::T reduceEdges(M map_f, R reduce_f) const {
    using T = typename R::T;
    auto D = parlay::delayed_seq<T>(n, [&](size_t i) {
      return get_vertex(i).out_neighbors().reduce(map_f, reduce_f);
    });
    return parlay::reduce(D, reduce_f);
  }

  // ======================= Constructors and fields  ========================
  symmetric_ptr_graph()
      : n(0),
        m(0),
        vertices(nullptr),
        edge_list_sizes(nullptr),
        vertex_weights(nullptr),
        deletion_fn([]() {}) {}

  symmetric_ptr_graph(size_t n, size_t m, vertex* _vertices,
                      std::function<void()> _deletion_fn,
                      vertex_weight_type* _vertex_weights = nullptr,
                      uintE* _edge_list_sizes = nullptr)
      : n(n),
        m(m),
        vertices(_vertices),
        edge_list_sizes(_edge_list_sizes),
        vertex_weights(_vertex_weights),
        deletion_fn(_deletion_fn) {}

  // Move constructor
  symmetric_ptr_graph(symmetric_ptr_graph&& other) noexcept {
    n = other.n;
    m = other.m;
    vertices = other.vertices;
    edge_list_sizes = other.edge_list_sizes;
    vertex_weights = other.vertex_weights;
    deletion_fn = std::move(other.deletion_fn);
    other.vertices = nullptr;
    other.edge_list_sizes = nullptr;
    other.vertex_weights = nullptr;
    other.deletion_fn = []() {};
  }

  // Move assignment
  symmetric_ptr_graph& operator=(symmetric_ptr_graph&& other) noexcept {
    n = other.n;
    m = other.m;
    vertices = other.vertices;
    edge_list_sizes = other.edge_list_sizes;
    vertex_weights = other.vertex_weights;
    deletion_fn();
    deletion_fn = std::move(other.deletion_fn);
    other.vertices = nullptr;
    other.edge_list_sizes = nullptr;
    other.vertex_weights = nullptr;
    other.deletion_fn = []() {};
  }

  // Copy constructor
  symmetric_ptr_graph(const symmetric_ptr_graph& other) {
    n = other.n;
    m = other.m;
    vertices = gbbs::new_array_no_init<vertex>(n);
    auto offsets = sequence<size_t>(n + 1);
    parallel_for(0, n, [&](size_t i) {
      vertices[i] = other.vertices[i];
      offsets[i] = vertices[i].out_degree();
    });
    offsets[n] = 0;
    size_t total_space = parlay::scan_inplace(make_slice(offsets));
    edge_type* E = gbbs::new_array_no_init<edge_type>(total_space);

    parallel_for(0, n, [&](size_t i) {
      size_t offset = offsets[i];
      auto map_f = [&](const uintE& u, const uintE& v, const W& wgh,
                       size_t ind) {
        E[offset + ind] = std::make_tuple(v, wgh);
      };
      // Copy neighbor data into E.
      vertices[i].out_neighbors().map_with_index(map_f);
      // Update this vertex's pointer to point to E.
      vertices[i].neighbors = E + offset;
    });

    deletion_fn = [this, E, total_space]() {
      gbbs::free_array(vertices, n);
      gbbs::free_array(E, total_space);
      if (vertex_weights != nullptr) {
        gbbs::free_array(vertex_weights, n);
      }
    };
    vertex_weights = nullptr;
    if (other.vertex_weights != nullptr) {
      vertex_weights = gbbs::new_array_no_init<vertex_weight_type>(n);
      parallel_for(
          0, n, [&](size_t i) { vertex_weights[i] = other.vertex_weights[i]; });
    }
  }

  ~symmetric_ptr_graph() { deletion_fn(); }

  // Note that observers recieve a handle to a vertex object which is only valid
  // so long as this graph's memory is valid.
  vertex get_vertex(uintE i) const { return vertices[i]; }

  // number of vertices in G
  size_t n;
  // number of edges in G
  size_t m;
  // pointer to array of vertex objects
  vertex* vertices;
  // pointer to array of vertex edge-list sizes---necessary if copying a
  // compressed graph in this representation.
  uintE* edge_list_sizes;
  // pointer to array of vertex weights
  vertex_weight_type* vertex_weights;

  // called to delete the graph
  std::function<void()> deletion_fn;
};

/* Compressed Sparse Row (CSR) based representation for asymmetric
 * graphs.  Note that the symmetric/asymmetric structures are pretty
 * similar, but defined separately. The purpose is to try and avoid
 * errors where an algorithm intended for symmetric graphs (e.g.,
 * biconnectivity) is not mistakenly called on a directed graph.
 *
 * Takes two template parameters:
 * 1) vertex_type: vertex template, parametrized by the weight type
 *    associated with each edge
 * 2) W: the edge weight template
 *
 * The graph is represented as an array of edges of type
 * vertex_type::edge_type, which is just a pair<uintE, W>.
 * */
template <template <class W> class vertex_type, class W>
struct asymmetric_graph {
  using vertex = vertex_type<W>;
  using weight_type = W;
  using edge_type = typename vertex::edge_type;
  using vertex_weight_type = double;

  // number of vertices in G
  size_t n;
  // number of edges in G
  size_t m;
  // called to delete the graph
  std::function<void()> deletion_fn;

  vertex_data* v_out_data;
  vertex_data* v_in_data;

  // Pointer to out-edges
  edge_type* out_edges;

  // Pointer to in-edges
  edge_type* in_edges;

  // Pointer to vertex weights
  vertex_weight_type* vertex_weights;

  vertex get_vertex(size_t i) const {
    return vertex(out_edges, v_out_data[i], in_edges, v_in_data[i], i);
  }

  asymmetric_graph()
      : n(0),
        m(0),
        deletion_fn([]() {}),
        v_out_data(nullptr),
        v_in_data(nullptr),
        out_edges(nullptr),
        in_edges(nullptr),
        vertex_weights(nullptr) {}

  asymmetric_graph(vertex_data* v_out_data, vertex_data* v_in_data, size_t n,
                   size_t m, std::function<void()> _deletion_fn,
                   edge_type* _out_edges, edge_type* _in_edges,
                   vertex_weight_type* _vertex_weights = nullptr)
      : n(n),
        m(m),
        deletion_fn(_deletion_fn),
        v_out_data(v_out_data),
        v_in_data(v_in_data),
        out_edges(_out_edges),
        in_edges(_in_edges),
        vertex_weights(_vertex_weights) {}

  // Move constructor
  asymmetric_graph(asymmetric_graph&& other) noexcept {
    n = other.n;
    m = other.m;
    v_out_data = other.v_out_data;
    v_in_data = other.v_in_data;
    out_edges = other.out_edges;
    in_edges = other.in_edges;
    vertex_weights = other.vertex_weights;
    deletion_fn = std::move(other.deletion_fn);
    other.v_out_data = nullptr;
    other.v_in_data = nullptr;
    other.out_edges = nullptr;
    other.in_edges = nullptr;
    other.vertex_weights = nullptr;
    other.deletion_fn = []() {};
  }

  // Move assignment
  asymmetric_graph& operator=(asymmetric_graph&& other) noexcept {
    n = other.n;
    m = other.m;
    v_out_data = other.v_out_data;
    v_in_data = other.v_in_data;
    out_edges = other.out_edges;
    in_edges = other.in_edges;
    vertex_weights = other.vertex_weights;
    deletion_fn();
    deletion_fn = std::move(other.deletion_fn);
    other.v_out_data = nullptr;
    other.v_in_data = nullptr;
    other.out_edges = nullptr;
    other.in_edges = nullptr;
    other.vertex_weights = nullptr;
    other.deletion_fn = []() {};
  }

  // Copy constructor
  asymmetric_graph(const asymmetric_graph& other) {
    debug(std::cout << "Copying asymmetric graph." << std::endl;);
    n = other.n;
    m = other.m;
    v_out_data = gbbs::new_array_no_init<vertex_data>(n);
    v_in_data = gbbs::new_array_no_init<vertex_data>(n);
    out_edges = gbbs::new_array_no_init<edge_type>(m);
    in_edges = gbbs::new_array_no_init<edge_type>(m);
    parallel_for(0, n, [&](size_t i) { v_out_data[i] = other.v_out_data[i]; });
    parallel_for(0, n, [&](size_t i) { v_in_data[i] = other.v_in_data[i]; });
    parallel_for(0, m, [&](size_t i) { out_edges[i] = other.out_edges[i]; });
    parallel_for(0, m, [&](size_t i) { in_edges[i] = other.in_edges[i]; });
    deletion_fn = [this]() {
      gbbs::free_array(v_out_data, n);
      gbbs::free_array(v_in_data, n);
      gbbs::free_array(out_edges, m);
      gbbs::free_array(in_edges, m);
      if (vertex_weights != nullptr) {
        gbbs::free_array(vertex_weights, n);
      }
    };
    vertex_weights = nullptr;
    if (other.vertex_weights != nullptr) {
      vertex_weights = gbbs::new_array_no_init<vertex_weight_type>(n);
      parallel_for(
          0, n, [&](size_t i) { vertex_weights[i] = other.vertex_weights[i]; });
    }
  }

  ~asymmetric_graph() { deletion_fn(); }

  template <class F>
  void mapEdges(F f, bool parallel_inner_map = true) const {
    parallel_for(0, n,
                 [&](size_t i) {
                   get_vertex(i).out_neighbors().map(f, parallel_inner_map);
                 },
                 1);
  }

  template <class F>
  void map_in_neighbors(size_t i, F f) const {
    auto f2 = [f](auto a, auto b, auto... args){return f(b,a,args...);};
    get_vertex(i).in_neighbors().map(f2, false);
  }
  template <class F>
  void map_out_neighbors(size_t i, F f) const {
    get_vertex(i).out_neighbors().map(f, false);
  }
  template <class F>
  void parallel_map_in_neighbors(size_t i, F f) const {
    auto f2 = [f](auto a, auto b, auto... args){return f(b,a,args...);};
    get_vertex(i).in_neighbors().map(f2, true);
  }
  template <class F>
  void parallel_map_out_neighbors(size_t i, F f) const {
    get_vertex(i).out_neighbors().map(f, true);
  }

  uintE out_degree(size_t i) const {
    return get_vertex(i).out_neighbors().degree;
  }

  uintE in_degree(size_t i) const {
    return get_vertex(i).in_neighbors().degree;
  }
};

// Similar to asymmetric_graph, but edges are not necessarily allocated
// consecutively.
template <template <class W> class vertex_type, class W>
struct asymmetric_ptr_graph {
  using vertex = vertex_type<W>;
  using weight_type = W;
  using edge_type = typename vertex::edge_type;
  using vertex_weight_type = double;

  // number of vertices in G
  size_t n;
  // number of edges in G
  size_t m;
  // pointer to array of vertex object
  vertex* vertices;
  // pointer to array of vertex weights
  vertex_weight_type* vertex_weights;

  // called to delete the graph
  std::function<void()> deletion_fn;

  vertex get_vertex(size_t i) const { return vertices[i]; }

  asymmetric_ptr_graph()
      : n(0),
        m(0),
        vertices(nullptr),
        deletion_fn([]() {}),
        vertex_weights(nullptr) {}

  asymmetric_ptr_graph(size_t n, size_t m, vertex* _vertices,
                       std::function<void()> _deletion_fn,
                       vertex_weight_type* _vertex_weights = nullptr)
      : n(n),
        m(m),
        vertices(_vertices),
        deletion_fn(_deletion_fn),
        vertex_weights(_vertex_weights) {}

  // Move constructor
  asymmetric_ptr_graph(asymmetric_ptr_graph&& other) noexcept {
    n = other.n;
    m = other.m;
    vertices = other.vertices;
    vertex_weights = other.vertex_weights;
    deletion_fn = std::move(other.deletion_fn);
    other.vertices = nullptr;
    other.vertex_weights = nullptr;
    other.deletion_fn = []() {};
  }

  // Move assignment
  asymmetric_ptr_graph& operator=(asymmetric_ptr_graph&& other) noexcept {
    n = other.n;
    m = other.m;
    vertices = other.vertices;
    vertex_weights = other.vertex_weights;
    deletion_fn();
    deletion_fn = std::move(other.deletion_fn);
    other.vertices = nullptr;
    other.vertex_weights = nullptr;
    other.deletion_fn = []() {};
  }

  // Copy constructor
  asymmetric_ptr_graph(const asymmetric_ptr_graph& other) {
    n = other.n;
    m = other.m;
    vertices = gbbs::new_array_no_init<vertex>(n);
    auto in_offsets = sequence<size_t>(n + 1);
    auto out_offsets = sequence<size_t>(n + 1);
    parallel_for(0, n, [&](size_t i) {
      vertices[i] = other.vertices[i];
      out_offsets[i] = vertices[i].out_degree();
      in_offsets[i] = vertices[i].in_degree();
    });
    in_offsets[n] = 0;
    out_offsets[n] = 0;

    size_t in_space = parlay::scan_inplace(make_slice(in_offsets));
    size_t out_space = parlay::scan_inplace(make_slice(out_offsets));
    edge_type* inE = gbbs::new_array_no_init<edge_type>(in_space);
    edge_type* outE = gbbs::new_array_no_init<edge_type>(out_space);

    parallel_for(0, n, [&](size_t i) {
      size_t out_offset = out_offsets[i];
      if (out_offsets[i + 1] != out_offset) {
        auto map_f = [&](const uintE& u, const uintE& v, const W& wgh,
                         size_t ind) {
          outE[out_offset + ind] = std::make_tuple(v, wgh);
        };
        // Copy neighbor data into E.
        vertices[i].out_neighbors().map_with_index(map_f);
        // Update this vertex's pointer to point to E.
        vertices[i].out_nghs = outE + out_offset;
      }

      size_t in_offset = in_offsets[i];
      if (in_offsets[i + 1] != in_offset) {
        auto map_f = [&](const uintE& u, const uintE& v, const W& wgh,
                         size_t ind) {
          inE[in_offset + ind] = std::make_tuple(v, wgh);
        };
        // Copy neighbor data into E.
        vertices[i].in_neighbors().map_with_index(map_f);
        // Update this vertex's pointer to point to E.
        vertices[i].in_nghs = inE + in_offset;
      }
    });

    deletion_fn = [this, inE, in_space, out_space, outE]() {
      gbbs::free_array(vertices, n);
      gbbs::free_array(inE, in_space);
      gbbs::free_array(outE, out_space);
      if (vertex_weights != nullptr) {
        gbbs::free_array(vertex_weights, n);
      }
    };
    vertex_weights = nullptr;
    if (other.vertex_weights != nullptr) {
      vertex_weights = gbbs::new_array_no_init<vertex_weight_type>(n);
      parallel_for(
          0, n, [&](size_t i) { vertex_weights[i] = other.vertex_weights[i]; });
    }
  }

  ~asymmetric_ptr_graph() { deletion_fn(); }

  template <class F>
  void mapEdges(F f, bool parallel_inner_map = true) const {
    parallel_for(0, n,
                 [&](size_t i) {
                   get_vertex(i).out_neighbors().map(f, parallel_inner_map);
                 },
                 1);
  }
};


template <class W, class set_type, bool inplace_ = false> struct set_wrapper {
  using weight_type = W;
  static constexpr bool binary = std::is_same_v<gbbs::empty, W>;
  static constexpr bool inplace = inplace_;
  static constexpr size_t inplace_cache_lines = 1;
  static constexpr size_t bytes_per_cache_line = 64;
  static_assert(
      inplace_cache_lines * bytes_per_cache_line > sizeof(set_type) || !inplace,
      "if we are storing edges in place we must have some bytes to do so");

  static constexpr int get_inplace_edge_count() {
    int bytes_we_have = inplace_cache_lines * bytes_per_cache_line;
    // to store the set
    bytes_we_have -= sizeof(set_type);
    // to store the count of elements in the in place array
    bytes_we_have -= 1;
    int elements_to_store;
    if constexpr (binary) {
      elements_to_store = bytes_we_have / sizeof(gbbs::uintE);
    } else {
      elements_to_store =
          bytes_we_have / sizeof(std::pair<gbbs::uintE, weight_type>);
    }

    return elements_to_store;
  }

  static constexpr int max_inplace_edge_count = get_inplace_edge_count();
  static_assert(max_inplace_edge_count < 256,
                "if we ever wanted to store more than 256 edges in place "
                "this would need to be fixed to not always be a uint8_t");

  static constexpr bool native_map = requires(const set_type &set) {
    set.map([](gbbs::uintE, weight_type) {});
  };

  static constexpr bool support_iteration = requires(const set_type &set) {
    ++set.begin();
    set.begin() != set.end();
  };

  static_assert(native_map || support_iteration);

  static constexpr bool native_map_early_exit = requires(const set_type &set) {
    set.map_early_exit([](gbbs::uintE, weight_type) {});
  };

  static constexpr bool native_parallel_map = requires(const set_type &set) {
    set.parallel_map([](gbbs::uintE, weight_type) {});
  };

  static constexpr bool native_parallel_map_early_exit =
      requires(const set_type &set) {
    set.parallel_map_early_exit([](gbbs::uintE, weight_type) {}, []() {});
  };

  static constexpr bool support_insert =
      requires(set_type &set, gbbs::uintE e) {
    set.insert(e);
    set.erase(e);
  };


  static constexpr bool native_insert_batch =
      requires(set_type &set, gbbs::uintE *es, size_t n) {
    set.insert_sorted_batch(es, n);
    set.remove_sorted_batch(es, n);
  };
  static constexpr bool insertable = support_insert | native_insert_batch;

  static void print_api() {
    std::cout << "### SET API ###\n";
    std::cout << "binary = " << binary << "\n";
    std::cout << "native map = " << native_map << "\n";
    std::cout << "iteration = " << support_iteration << "\n";
    std::cout << "native map_early_exit = " << native_map_early_exit << "\n";
    std::cout << "native parallel_map = " << native_parallel_map << "\n";
    std::cout << "native native_parallel_map_early_exit = "
              << native_parallel_map_early_exit << "\n";
    if constexpr (inplace) {
      std::cout << "storing up to " << max_inplace_edge_count
                << " edges in place\n";
    }
    if constexpr (insertable) {
      std::cout << "supports modification\n";
      if constexpr (native_insert_batch) {
        std::cout << "supports insert_batch\n";
      }
    } else {
      std::cout << "static structure cannot be updated\n";
    }
  }

  size_t size() const {
    if constexpr (inplace) {
      return inplace_count + set.size();
    } else {
      return set.size();
    }
  }

  // template <class InputIt1, class InputIt2>
  // std::array<inplace_array_element_type, max_inplace_edge_count>
  // merge(InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2) {
  //   for (; first1 != last1; ++d_first) {
  //     if (first2 == last2)
  //       return std::copy(first1, last1, d_first);

  //     if (*first2 < *first1)
  //       *d_first = *first2++;
  //     else {
  //       *d_first = *first1;
  //       if (!(*first1 < *first2))
  //         ++first2;
  //       ++first1;
  //     }
  //   }
  //   return std::copy(first2, last2, d_first);
  // }

  void sorted_batch_insert_no_inplace(gbbs::uintE *es, size_t n) {
    if constexpr (native_insert_batch) {
      set.insert_sorted_batch(es, n);
    } else if constexpr (support_insert) {
      for (size_t i = 0; i < n; i++) {
        set.insert(es[i]);
      }
    }
  }

   void insert_sorted_batch(gbbs::uintE *es, size_t n) {
    // std::cout << "batch\n";
    // for (int i = 0; i < n; i++) {
    //   std::cout << es[i] << ", ";
    // }
    // std::cout << "\n";
    // std::cout << "set\n";
    // map<false>([](auto elem, weight_type w) {
    //   std::cout << elem << ", ";
    // });
    // std::cout << "\n";
    if constexpr (inplace) {
      std::array<inplace_array_element_type, max_inplace_edge_count>
          new_inplace = {};
      auto inplace_it = inplace_array.data();
      auto inplace_end = inplace_array.data() + inplace_count;
      gbbs::uintE *batch_end = es + n;
      int new_inplace_count = 0;
      bool done_in_place = false;

      while (inplace_it != inplace_end && es != batch_end && new_inplace_count < max_inplace_edge_count) {
        if (*inplace_it < *es) {
          new_inplace[new_inplace_count] = *inplace_it;
          inplace_it++;
          new_inplace_count++;
          
        } else if (*es < *inplace_it) {
          new_inplace[new_inplace_count] = *es;
          es++;
          new_inplace_count++;
          
        } else {
          // they are equal
          new_inplace[new_inplace_count] = *inplace_it;
          inplace_it++;
          es++;
          new_inplace_count++;
        }
      }

      // there is still data in one of them that needs to fill the inplace
      if (new_inplace_count < max_inplace_edge_count) {
        if (inplace_it < inplace_end) {
          while (new_inplace_count < max_inplace_edge_count && inplace_it < inplace_end) {
            new_inplace[new_inplace_count] = *inplace_it;
            new_inplace_count++;
            inplace_it++;
          }
          if (inplace_it != inplace_end) {
            // may still be some data in the old inplace
            sorted_batch_insert_no_inplace(inplace_it, inplace_end - inplace_it);
          }
          inplace_array = new_inplace;
          inplace_count = new_inplace_count;
          return;
        } else {
          while (new_inplace_count < max_inplace_edge_count && es < batch_end) {
            new_inplace[new_inplace_count] = *es;
            new_inplace_count++;
            es++;
          }
          if (es != batch_end) {
            // may still be some data in the batch
            sorted_batch_insert_no_inplace(es, batch_end - es);
          }
          inplace_array = new_inplace;
          inplace_count = new_inplace_count;
          return;
        }
      }
      // we are done with the inplace, there may be data in both

      if (inplace_it < inplace_end) {
        //if something got kicked out of the inplace, there must be room for it at the front of the batch
        size_t num_kicked = inplace_end - inplace_it;
        sorted_batch_insert_no_inplace(inplace_it, num_kicked);
      }
      inplace_array = new_inplace;
      inplace_count = new_inplace_count;
      n = batch_end - es;
      if (n == 0) {
        return;
      }
    }
    sorted_batch_insert_no_inplace(es, n);
  }

  void sorted_batch_remove_no_inplace(gbbs::uintE *es, size_t n) {
    if constexpr (native_insert_batch) {
      set.remove_sorted_batch(es, n);
    } else if constexpr (support_insert) {
      for (size_t i = 0; i < n; i++) {
        set.erase(es[i]);
      }
    }
  }

  void remove_sorted_batch(gbbs::uintE *es, size_t n) {

    if constexpr (inplace) {
      size_t num_written = 0;
      auto batch_end = es+n;
      for (int i = 0; i < inplace_count; i++) {
        auto it = std::lower_bound(es, batch_end, inplace_array[i]);
        if (it != batch_end && *it != inplace_array[i]) {
          inplace_array[num_written] = inplace_array[i];
          num_written++;
        }
      }
      for (int i = num_written; i < inplace_count; i++) {
        inplace_array[i] = {};
      }
      inplace_count = num_written;
      // now we need to refil the inplace with data from the set
      // this is only correct for sorted sets, but there are lots of other issues that arise from unsorted structures with the inplace during inserts and deletes, so we ignore them for now
      if (inplace_count < max_inplace_edge_count) {
        map_early_exit<true>([&](auto elem, [[maybe_unused]] weight_type w) {
          if (inplace_count < max_inplace_edge_count) {
              inplace_array[inplace_count] = elem;
              inplace_count++;
          }
          return inplace_count == max_inplace_edge_count;
        });
        int num_taken = inplace_count - num_written;
        sorted_batch_remove_no_inplace(inplace_array.data() + num_written, num_taken);
      }
      n = batch_end - es;
    }
    sorted_batch_remove_no_inplace(es, n);
  }

  template <bool inplace_done = false, class F> void map(F f) const {
    if constexpr (inplace && !inplace_done) {
      for (int i = 0; i < inplace_count; i++) {
        if constexpr (binary) {
          weight_type w{};
          f(inplace_array[i], w);
        } else {
          f(inplace_array[i].first, inplace_array[i].second);
        }
      }
    }
    if constexpr (native_map) {
      set.map(f);
    } else if constexpr (support_iteration) {
      if constexpr (binary) {
        for (const auto &other : set) {
          weight_type w{};
          f(other, w);
        }
      } else {
        for (const auto &[other, weight] : set) {
          f(other, weight);
        }
      }
    }
  }

  template <bool inplace_done = false, class F> void map_early_exit(F f) const {
    if constexpr (inplace && !inplace_done) {
      for (int i = 0; i < inplace_count; i++) {
        if constexpr (binary) {
          weight_type w{};
          if (f(inplace_array[i], w)) {
            return;
          }
        } else {
          if (f(inplace_array[i].first, inplace_array[i].second)) {
            return;
          }
        }
      }
    }
    if constexpr (native_map_early_exit) {
      set.map_early_exit(f);
    } else if constexpr (support_iteration) {
      if constexpr (binary) {
        for (const auto &other : set) {
          weight_type w{};
          if (f(other, w)) {
            break;
          }
        }
      } else {
        for (const auto &[other, weight] : set) {
          if (f(other, weight)) {
            break;
          }
        }
      }
    } else {
      map<true>(f);
    }
  }

  template <bool inplace_done = false, class F> void parallel_map(F f) const {
    if constexpr (inplace && !inplace_done) {
      for (int i = 0; i < inplace_count; i++) {
        if constexpr (binary) {
          weight_type w{};
          f(inplace_array[i], w);
        } else {
          f(inplace_array[i].first, inplace_array[i].second);
        }
      }
    }
    if constexpr (native_parallel_map) {
      set.parallel_map(f);
    } else {
      map<true>(f);
    }
  }

  template <class F, typename Block_F = std::nullptr_t>
  void parallel_map_early_exit(F f, Block_F block_check = {}) const {
    if constexpr (inplace) {
      for (int i = 0; i < inplace_count; i++) {
        if constexpr (binary) {
          weight_type w{};
          if (f(inplace_array[i], w)) {
            return;
          }
        } else {
          if (f(inplace_array[i].first, inplace_array[i].second)) {
            return;
          }
        }
      }
    }
    if constexpr (native_parallel_map_early_exit) {
      set.parallel_map_early_exit(f, block_check);
    } else if constexpr (native_map_early_exit) {
      map_early_exit<true>(f);
    } else if constexpr (native_parallel_map) {
      parallel_map<true>(f);
    } else {
      map<true>(f);
    }
  }

  set_wrapper() = default;

  set_wrapper(auto begin, auto end) {
    // TODO(wheatman) handle other types of creation
    if constexpr (inplace) {
      for (int i = 0; i < max_inplace_edge_count; i++) {
        if (begin < end) {
          if constexpr (binary) {
            inplace_array[i] = std::get<0>(*begin);
          } else {
            inplace_array[i] = {std::get<0>(*begin), std::get<1>(*begin)};
          }
          ++begin;
          inplace_count = i + 1;
        }
      }
    }

    // this is done so we can call the constructor after we deal with the
    // possible inplace above. This assume that we have an efficient move
    // operator The other set should still be empty, so even buggy move
    // operators should still work fine
    // however we need to convert the range from tuples to elements if we are unweighted
    auto range =
        std::views::transform(std::ranges::subrange(begin, end), [](auto elem) {
          if constexpr (binary) {
            return std::get<0>(elem);
          } else {
            return elem;
          }
        });
    set_type set2(range.begin(), range.end());
    set = std::move(set2);
  }

  // Move constructor
  set_wrapper(set_wrapper &&other) noexcept {
    if constexpr (inplace) {
      inplace_array = other.inplace_array;
      inplace_count = other.inplace_count;
    }
    set = other.set;
  }

  // Copy constructor
  set_wrapper(const set_wrapper &other) noexcept {
    if constexpr (inplace) {
      inplace_array = other.inplace_array;
      inplace_count = other.inplace_count;
    }
    set = other.set;
  }

  // Move assignment
  set_wrapper &operator=(set_wrapper &&other) noexcept {
    if constexpr (inplace) {
      inplace_array = other.inplace_array;
      inplace_count = other.inplace_count;
    }
    set = other.set;
    return *this;
  }

  // copy assignment
  set_wrapper &operator=(const set_wrapper &other) noexcept {
    if constexpr (inplace) {
      inplace_array = other.inplace_array;
      inplace_count = other.inplace_count;
    }
    set = other.set;
    return *this;
  }

  ~set_wrapper() = default;

private:
  using inplace_array_element_type =
      std::conditional<binary, gbbs::uintE,
                       std::pair<gbbs::uintE, weight_type>>::type;
  set_type set = {};
  [[no_unique_address]] typename std::conditional<
      inplace, std::array<inplace_array_element_type, max_inplace_edge_count>,
      gbbs::empty>::type inplace_array;

  [[no_unique_address]]
  typename std::conditional<inplace, int, gbbs::empty>::type inplace_count = {};
};


//  Adjencency List based representation for symmetric graphs.
//  Takes tree template parameters:
//  1) vertex_type: vertex template, parametrized by the weight type associated
//  with each edge
//  2) W: the edge weight template
//  3) set: the set implementation to use, will be wrapper in a set_wrapper
template <template <class W> class vertex_type, class W, class set, bool inplace = false>
struct symmetric_set_graph {
  using vertex = vertex_type<W>;
  using weight_type = W;
  static constexpr bool binary = std::is_same_v<gbbs::empty, W>;
  using edge_type = typename vertex::edge_type;
  using vertex_weight_type = double;
  using set_wrapper_type = set_wrapper<W, set, inplace>;
  static constexpr bool support_inserts = set_wrapper_type::insertable;

  size_t N() const { return nodes.size(); }

  uintE degree(size_t i) const {
    return nodes[i].size();
  }


  template <class F>
  void map_neighbors(size_t i, F f) const {
    nodes[i].map([&](gbbs::uintE other, weight_type w) {
      return f(i, other, w);
    });
  }

  template <class F>
  void map_neighbors_early_exit(size_t i, F f) const {
    nodes[i].map_early_exit([&](gbbs::uintE other, weight_type w) {
      return f(i, other, w);
    });
  }

  template <class F>
  void parallel_map_neighbors(size_t i, F f) const {
    nodes[i].parallel_map([&](gbbs::uintE other, weight_type w) {
      return f(i, other, w);
    });
  }

  template <class F, typename Block_F = std::nullptr_t>
  void parallel_map_neighbors_early_exit(size_t i, F f, Block_F block_check = {}) const {
    nodes[i].parallel_map_early_exit([&](gbbs::uintE other, weight_type w) {
      return f(i, other, w);
    }, block_check);
  }


  void insert_sorted_grouped_batch(parlay::sequence<std::pair<gbbs::uintE, parlay::sequence<gbbs::uintE>>> & els) {
    parlay::parallel_for(0, els.size(), [&](size_t i) {
      if (!els[i].second.empty()) {
        nodes[els[i].first].insert_sorted_batch(els[i].second.data(),
                                          els[i].second.size());
      }
    });
  }

  void remove_sorted_grouped_batch(parlay::sequence<std::pair<gbbs::uintE, parlay::sequence<gbbs::uintE>>> & els) {
    parlay::parallel_for(0, els.size(), [&](size_t i) {
      if (!els[i].second.empty()) {
        nodes[els[i].first].remove_sorted_batch(els[i].second.data(),
                                          els[i].second.size());
      }
    });
  }

  // ======================= Constructors and fields  ========================
  symmetric_set_graph()
      :vertex_weights(nullptr) {}

  // for now use the same constructor as the exsisting one for ease, but probably write a more general constuctor for everyone to use later
  // TODO(wheatman) figure out what to do with _deletion_fn, probably should just use unique pointer
  symmetric_set_graph(vertex_data* v_data, size_t n, size_t m,
                  std::function<void()> _deletion_fn, edge_type* _e0,
                  vertex_weight_type* _vertex_weights = nullptr)
      : nodes(n), vertex_weights(_vertex_weights), deletion_fn(_deletion_fn)
        {
          //TODO(wheatman) this is just for debugging help to check that the correct api is generated for different sets)
          set_wrapper_type::print_api();
          parallel_for(0, n, [&](size_t i){
            nodes[i] = set_wrapper_type(&_e0[v_data[i].offset], &_e0[v_data[i].offset+v_data[i].degree]);
          });
        }

  // Move constructor
  symmetric_set_graph(symmetric_set_graph&& other) noexcept {
    nodes = other.nodes;
    vertex_weights = other.vertex_weights;
    other.vertex_weights = nullptr;
    deletion_fn = std::move(other.deletion_fn);
  }

  // Move assignment
  symmetric_set_graph& operator=(symmetric_set_graph&& other) noexcept {
    nodes = other.nodes;
    vertex_weights = other.vertex_weights;
    other.vertex_weights = nullptr;
    deletion_fn = std::move(other.deletion_fn);
  }

  // Copy constructor
  symmetric_set_graph(const symmetric_set_graph& other) {
    debug(std::cout << "Copying symmetric graph." << std::endl;);
    nodes = other.nodes;
    vertex_weights = nullptr;
    if (other.vertex_weights != nullptr) {
      vertex_weights = gbbs::new_array_no_init<vertex_weight_type>(nodes.size());
      parallel_for(
          0, nodes.size(), [&](size_t i) { vertex_weights[i] = other.vertex_weights[i]; });
    }
    deletion_fn = [&]() {
      if (vertex_weights != nullptr) {
        gbbs::free_array(vertex_weights, nodes.size());
      }
    };
  }

  ~symmetric_set_graph() {
      deletion_fn();
   }

  // Graph Data
  std::vector<set_wrapper_type> nodes;
  vertex_weight_type* vertex_weights;
  // called to delete the graph
  std::function<void()> deletion_fn;

};


// Mutates (sorts) the underlying array A containing a black-box description of
// an edge of typename A::value_type. The caller provides functions GetU, GetV,
// and GetW which extract the u, v, and weight of a (u,v,w) edge respective (if
// the edge is a std::tuple<uinte, uintE, W> this is just get<0>, ..<1>, ..<2>
// respectively.
// e.g.:
//   using edge = std::tuple<uintE, uintE, W>;
//   auto get_u = [&] (const edge& e) { return std::get<0>(e); };
//   auto get_v = [&] (const edge& e) { return std::get<1>(e); };
//   auto get_w = [&] (const edge& e) { return std::get<2>(e); };
//   auto G = sym_graph_from_edges<W>(coo1, get_u, get_v, get_w, 10, false);
template <class graph, class EdgeSeq, class GetU, class GetV, class GetW>
static inline graph sym_graph_from_edges(
    EdgeSeq& A, size_t n, GetU&& get_u, GetV&& get_v, GetW&& get_w,
    bool is_sorted = false) {
  using Wgh = typename graph::weight_type;
  using vertex = symmetric_vertex<Wgh>;
  using edge_type = typename vertex::edge_type;
  size_t m = A.size();

  if (m == 0) {
    if (n == 0) {
      std::function<void()> del = []() {};
      return graph(nullptr, 0, 0, del,
                                                    nullptr);
    } else {
      auto v_data = gbbs::new_array_no_init<vertex_data>(n);
      parallel_for(0, n, [&](size_t i) {
        v_data[i].offset = 0;
        v_data[i].degree = 0;
      });
      return graph(
          v_data, n, 0, [=]() { gbbs::free_array(v_data, n); }, nullptr);
    }
  }

  if (!is_sorted) {
    parlay::integer_sort_inplace(make_slice(A), get_u);
  }

  auto starts = sequence<uintT>(n + 1, (uintT)0);

  using neighbor = std::tuple<uintE, Wgh>;
  auto edges = gbbs::new_array_no_init<neighbor>(m);
  parallel_for(0, m, [&](size_t i) {
    if (i == 0 || (get_u(A[i]) != get_u(A[i - 1]))) {
      starts[get_u(A[i])] = i;
    }
    if (i != (m - 1)) {
      uintE our_vtx = get_u(A[i]);
      uintE next_vtx = get_u(A[i + 1]);
      if (our_vtx != next_vtx && (our_vtx + 1 != next_vtx)) {
        parallel_for(our_vtx + 1, next_vtx,
                     [&](size_t k) { starts[k] = i + 1; });
      }
    }
    if (i == (m - 1)) { /* last edge */
      parallel_for(get_u(A[i]) + 1, starts.size(),
                   [&](size_t j) { starts[j] = m; });
    }
    edges[i] = std::make_tuple(get_v(A[i]), get_w(A[i]));
  });

  auto v_data = gbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&](size_t i) {
    uintT o = starts[i];
    v_data[i].offset = o;
    v_data[i].degree = (uintE)(((i == (n - 1)) ? m : starts[i + 1]) - o);
  });
  return graph(
      v_data, n, m,
      [=]() {
        gbbs::free_array(v_data, n);
        gbbs::free_array(edges, m);
      },
      (edge_type *)edges);
}

template <class graph>
static inline graph sym_graph_from_edges(
    sequence<std::tuple<uintE, uintE, typename graph::weight_type>>& A, size_t n,
    bool is_sorted = false) {
  using edge = std::tuple<uintE, uintE, typename graph::weight_type>;
  auto get_u = [&](const edge& e) { return std::get<0>(e); };
  auto get_v = [&](const edge& e) { return std::get<1>(e); };
  auto get_w = [&](const edge& e) { return std::get<2>(e); };
  return graph_implementations::template sym_graph_from_edges<graph>(A, n, get_u, get_v, get_w, is_sorted);
}

template <class Wgh, class EdgeSeq, class GetU, class GetV, class GetW>
std::tuple<uintE, Wgh>* get_edges(EdgeSeq& A, sequence<uintT>& starts, size_t m,
                                  const GetU& get_u, const GetV& get_v,
                                  const GetW& get_w) {
  using neighbor = std::tuple<uintE, Wgh>;
  auto edges = gbbs::new_array_no_init<neighbor>(m);
  parallel_for(0, m, [&](size_t i) {
    if (i == 0 || (get_u(A[i]) != get_u(A[i - 1]))) {
      starts[get_u(A[i])] = i;
    }
    if (i != (m - 1)) {
      uintE our_vtx = get_u(A[i]);
      uintE next_vtx = get_u(A[i + 1]);
      if (our_vtx != next_vtx && (our_vtx + 1 != next_vtx)) {
        parallel_for(our_vtx + 1, next_vtx,
                     [&](size_t k) { starts[k] = i + 1; });
      }
    }
    if (i == (m - 1)) { /* last edge */
      parallel_for(get_u(A[i]) + 1, starts.size(),
                   [&](size_t j) { starts[j] = m; });
    }
    edges[i] = std::make_tuple(get_v(A[i]), get_w(A[i]));
  });
  return edges;
}

// Mutates (sorts) the underlying array A containing a black-box description of
// an edge of typename A::value_type. The caller provides functions GetU, GetV,
// and GetW which extract the u, v, and weight of a (u,v,w) edge respective (if
// the edge is a std::tuple<uinte, uintE, W> this is just get<0>, ..<1>, ..<2>
// respectively.
// e.g.:
//   using edge = std::tuple<uintE, uintE, W>;
//   auto get_u = [&] (const edge& e) { return std::get<0>(e); };
//   auto get_v = [&] (const edge& e) { return std::get<1>(e); };
//   auto get_w = [&] (const edge& e) { return std::get<2>(e); };
//   auto G = asym_graph_from_edges<W>(coo1, get_u, get_v, get_w, 10, false);
template <class Wgh, class EdgeSeq, class GetU, class GetV, class GetW>
static inline asymmetric_graph<asymmetric_vertex, Wgh> asym_graph_from_edges(
    EdgeSeq& A, size_t n, GetU&& get_u, GetV&& get_v, GetW&& get_w,
    bool is_sorted = false) {
  using vertex = asymmetric_vertex<Wgh>;
  using edge_type = typename vertex::edge_type;
  size_t m = A.size();

  if (m == 0) {
    if (n == 0) {
      std::function<void()> del = []() {};
      return asymmetric_graph<asymmetric_vertex, Wgh>(nullptr, nullptr, 0, 0,
                                                      del, nullptr, nullptr);
    } else {
      auto v_in_data = gbbs::new_array_no_init<vertex_data>(n);
      auto v_out_data = gbbs::new_array_no_init<vertex_data>(n);
      parallel_for(0, n, [&](size_t i) {
        v_in_data[i].offset = 0;
        v_in_data[i].degree = 0;

        v_out_data[i].offset = 0;
        v_out_data[i].degree = 0;
      });
      return asymmetric_graph<asymmetric_vertex, Wgh>(
          v_out_data, v_in_data, n, 0,
          [=]() {
            gbbs::free_array(v_out_data, n);
            gbbs::free_array(v_in_data, n);
          },
          nullptr, nullptr);
    }
  }

  // flip to create the in-edges
  auto I = sequence<typename EdgeSeq::value_type>::from_function(
      A.size(), [&](size_t i) {
        using T = typename EdgeSeq::value_type;
        auto e = A[i];
        return T(get_v(e), get_u(e), get_w(e));
      });

  if (!is_sorted) {
    parlay::integer_sort_inplace(make_slice(A), get_u);
    parlay::integer_sort_inplace(make_slice(I), get_u);
  }

  auto in_starts = sequence<uintT>(n + 1, (uintT)0);
  auto out_starts = sequence<uintT>(n + 1, (uintT)0);

  auto in_edges = get_edges<Wgh>(I, in_starts, m, get_u, get_v, get_w);
  auto out_edges = get_edges<Wgh>(A, out_starts, m, get_u, get_v, get_w);

  auto in_v_data = gbbs::new_array_no_init<vertex_data>(n);
  auto out_v_data = gbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n,
               [&](size_t i) {
                 uintT in_o = in_starts[i];
                 in_v_data[i].offset = in_o;
                 in_v_data[i].degree =
                     (uintE)(((i == (n - 1)) ? m : in_starts[i + 1]) - in_o);

                 uintT out_o = out_starts[i];
                 out_v_data[i].offset = out_o;
                 out_v_data[i].degree =
                     (uintE)(((i == (n - 1)) ? m : out_starts[i + 1]) - out_o);
               },
               1024);
  return asymmetric_graph<asymmetric_vertex, Wgh>(
      out_v_data, in_v_data, n, m,
      [=]() {
        gbbs::free_array(in_v_data, n);
        gbbs::free_array(out_v_data, n);
        gbbs::free_array(in_edges, m);
        gbbs::free_array(out_edges, m);
      },
      (edge_type*)out_edges, (edge_type*)in_edges);
}

template <class Wgh>
static inline asymmetric_graph<asymmetric_vertex, Wgh> asym_graph_from_edges(
    sequence<std::tuple<uintE, uintE, Wgh>>& A, size_t n,
    bool is_sorted = false) {
  using edge = std::tuple<uintE, uintE, Wgh>;
  auto get_u = [&](const edge& e) { return std::get<0>(e); };
  auto get_v = [&](const edge& e) { return std::get<1>(e); };
  auto get_w = [&](const edge& e) { return std::get<2>(e); };
  return asym_graph_from_edges<Wgh>(A, n, get_u, get_v, get_w, is_sorted);
}

} // namespace graph_implementations

template <bool degree_, bool M_, bool parallel_map_, bool map_early_exit_, bool parallel_map_early_exit_, bool map_range_, bool store_M_, bool store_degrees_>
class Graph_API_use {
  public:
  static constexpr bool degree = degree_;
  static constexpr bool M = M_;
  static constexpr bool parallel_map = parallel_map_;
  static constexpr bool map_early_exit = map_early_exit_;
  static constexpr bool parallel_map_early_exit = parallel_map_early_exit_;
  static constexpr bool store_M = store_M_;
  static constexpr bool store_degrees = store_degrees_;
  static constexpr bool map_range = map_range_;
};

using full_api = Graph_API_use<true, true, true, true, false, true, true, true>;
using minimal_api = Graph_API_use<false, false, false, false, false, false, false, false>;
using no_early_exit = Graph_API_use<true, true, true, false, false, true, true, true>;
using no_parallel_map = Graph_API_use<true, true, false, true, false, true, true, true>;
using no_degree = Graph_API_use<false, true, true, true, true, true, true, false>;
using no_M = Graph_API_use<true, false, true, true, true, true, true, false>;

template <class Representation, bool symmetric, class API = full_api> class Graph {
  Representation g;
  using api = API;

  static constexpr bool support_native_degree =
      requires(const Representation &g, uint64_t i) {
    g.degree(i);
  };

  static constexpr bool support_native_out_degree =
      requires(const Representation &g, uint64_t i) {
    g.out_degree(i);
  };
  static constexpr bool support_native_in_degree =
      requires(const Representation &g, uint64_t i) {
    g.in_degree(i);
  };
  static constexpr bool support_native_degree_directred =
      support_native_out_degree && support_native_in_degree;

  static constexpr bool n_member = requires(const Representation &g) { g.n; };
  static constexpr bool n_function = requires(const Representation &g) { g.N(); };
  static_assert(n_member || n_function,
                "the underlying graph representation must support getting "
                "the number of nodes\n");

  static constexpr bool example_map_func(uintE a, uintE b, Representation::weight_type w) {
    return true;
  }
  static constexpr bool null_func() {return true;}

  static constexpr bool support_map_neighbor =
      requires(const Representation &g, size_t i) {
    g.map_neighbors(i, example_map_func);
  };
  

  static constexpr bool support_map_in_neighbor =
      requires(const Representation &g, size_t i) {
    g.map_in_neighbors(i, example_map_func);
  };
  

  static constexpr bool support_map_out_neighbor =
      requires(const Representation &g, size_t i) {
    g.map_out_neighbors(i, example_map_func);
  };
  

  static constexpr bool support_parallel_map_neighbors =
      requires(const Representation &g, size_t i) {
    g.parallel_map_neighbors(i, example_map_func);
  };
  static constexpr bool support_parallel_map_in_neighbors =
      requires(const Representation &g, size_t i) {
    g.parallel_map_in_neighbors(i, example_map_func);
  };
  static constexpr bool support_parallel_map_out_neighbors =
      requires(const Representation &g, size_t i) {
    g.parallel_map_out_neighbors(i, example_map_func);
  };
  static constexpr bool support_map_neighbors_early_exit =
      requires(const Representation &g, size_t i) {
    g.map_neighbors_early_exit(i, example_map_func);
  };
  static constexpr bool support_map_in_neighbors_early_exit =
      requires(const Representation &g, size_t i) {
    g.map_in_neighbors_early_exit(i, example_map_func);
  };
  static constexpr bool support_map_out_neighbors_early_exit =
      requires(const Representation &g, size_t i) {
    g.map_out_neighbors_early_exit(i, example_map_func);
  };
  static constexpr bool support_parallel_map_neighbors_early_exit1 =
      requires(const Representation &g, size_t i) {
    g.parallel_map_neighbors_early_exit(i, example_map_func);
  };
  static constexpr bool support_parallel_map_neighbors_early_exit2 =
      requires(const Representation &g, size_t i) {
    g.parallel_map_neighbors_early_exit(i, example_map_func, null_func);
  };
  static constexpr bool support_parallel_map_neighbors_early_exit =
      support_parallel_map_neighbors_early_exit1 ||
      support_parallel_map_neighbors_early_exit2;

  static constexpr bool support_parallel_map_in_neighbors_early_exit1 =
      requires(const Representation &g, size_t i) {
    g.parallel_map_in_neighbors_early_exit(i, example_map_func);
  };
  static constexpr bool support_parallel_map_in_neighbors_early_exit2 =
      requires(const Representation &g, size_t i) {
    g.parallel_map_in_neighbors_early_exit(i, example_map_func, null_func);
  };
  static constexpr bool support_parallel_map_in_neighbors_early_exit =
      support_parallel_map_in_neighbors_early_exit1 ||
      support_parallel_map_in_neighbors_early_exit2;
  static constexpr bool support_parallel_map_out_neighbors_early_exit1 =
      requires(const Representation &g, size_t i) {
    g.parallel_map_out_neighbors_early_exit(i, example_map_func);
  };
  static constexpr bool support_parallel_map_out_neighbors_early_exit2 =
      requires(const Representation &g, size_t i) {
    g.parallel_map_out_neighbors_early_exit(i, example_map_func, null_func);
  };
  static constexpr bool support_parallel_map_out_neighbors_early_exit =
      support_parallel_map_out_neighbors_early_exit1 ||
      support_parallel_map_out_neighbors_early_exit2;
  static constexpr bool m_member = requires(const Representation &g) { g.m; };
  static constexpr bool m_function = requires(const Representation &g) {
    g.M();
  };
  static constexpr bool storing_m =
      (!m_function) && (!m_member) && api::store_M && api::M;

  static constexpr bool storing_degrees = (!support_native_degree) &&
                                          api::store_degrees && symmetric &&
                                          api::degree;

  static constexpr bool storing_directed_degrees =
      (!support_native_degree_directred) && api::store_degrees && !symmetric &&
      api::degree;



  static constexpr bool set_wrapper_supports_inserts_() {
    constexpr bool support_insert_member = requires(Representation & g) {
      Representation::insertable;
    };
    if constexpr(support_insert_member) {
      return Representation::insertable;
    } else {
      return true;
    }
  }


  std::conditional<storing_m, size_t, gbbs::empty>::type stored_m;
  std::conditional<storing_degrees, parlay::sequence<uintE>, gbbs::empty>::type
      stored_degrees;
  std::conditional<storing_directed_degrees, parlay::sequence<uintE>,
                   gbbs::empty>::type stored_in_degrees;
  std::conditional<storing_directed_degrees, parlay::sequence<uintE>,
                   gbbs::empty>::type stored_out_degrees;

  static constexpr bool needs_finalize_() {
    if constexpr (storing_m) {
      return true;
    }
    if constexpr (symmetric) {
      if constexpr (storing_degrees) {
        return true;
      }
    } else {
      if constexpr (storing_directed_degrees) {
        return true;
      }
    }
    return false;
  }


  uintE degree_no_store(size_t i) const {
    static_assert(symmetric);
    if constexpr (support_native_degree && api::degree) {
      return g.degree(i);
    } else {
      // TODO(wheatman) if we know that this is serial than we don't need a
      // reducer here
      Reducer_sum<size_t> count;
      g.map_neighbors(i,
                      [&count]([[maybe_unused]] auto... args) { count += 1; });
      return count;
    }
  }

    uintE out_degree_no_store(size_t i) const {
    if constexpr (support_native_degree_directred && api::degree) {
      return g.out_degree(i);
    } else {
      if constexpr (symmetric) {
        return degree_no_store(i);
      } else {
        // TODO(wheatman) if we know that this is serial than we don't need a
        // reducer here
        Reducer_sum<size_t> count;
        g.map_out_neighbors(
            i, [&count]([[maybe_unused]] auto... args) { count += 1; });
        return count;
      }
    }
  }

  uintE in_degree_no_store(size_t i) const {
    constexpr bool support_native =
        requires(const Representation &g, uint64_t i) {
      g.in_degree(i);
    };
    if constexpr (support_native && api::degree) {
      return g.in_degree(i);
    } else {
      if constexpr (symmetric) {
        return degree_no_store(i);
      } else {
        // TODO(wheatman) if we know that this is serial than we don't need a
        // reducer here
        Reducer_sum<size_t> count;
        g.map_in_neighbors(
            i, [&count]([[maybe_unused]] auto... args) { count += 1; });
        return count;
      }
    }
  }

public:
  static constexpr bool support_insert_sorted_batch = requires(
      Representation & g, std::tuple<gbbs::uintE, gbbs::uintE> *es, size_t n) {
    g.insert_sorted_batch(es, n);
    g.remove_sorted_batch(es, n);
  };


  static constexpr bool support_insert_grouped_batch = requires(
      Representation & g, parlay::sequence<std::pair<gbbs::uintE, parlay::sequence<gbbs::uintE>>> & els) {
    g.insert_sorted_grouped_batch(els);
    g.remove_sorted_grouped_batch(els);
  };
  static constexpr bool support_insert_batch =
      (support_insert_sorted_batch || support_insert_grouped_batch) &&
      set_wrapper_supports_inserts_();
  static constexpr bool needs_finalize = needs_finalize_();
  static void print_api() {
    std::cout << "### GRAPH API ###\n";
    if constexpr (n_member) {
      std::cout << "n is suppoted as a member\n";
    } else {
      static_assert(n_function);
      std::cout << "n is suppoted as a function\n";
    }
    if constexpr (m_member) {
      std::cout << "m is suppoted as a member\n";
    } else if constexpr (m_function) {
      std::cout << "m is suppoted as a function\n";
    } else if constexpr (storing_m) {
      std::cout << "m is stored\n";
    } else {
      std::cerr << "m will be calculated as needed\nTHIS IS VERY INEFFICIENT\n";
    }
    if constexpr (symmetric) {
      std::cout << "the graph is symmetric\n";
      if constexpr (support_native_degree) {
        std::cout << "degree is supported\n";
      } else if constexpr (storing_degrees) {
        std::cout << "degrees are stored\n";
      } else {
        std::cerr << "degrees will be calculated as needed\nTHIS IS VERY "
                     "INEFFICIENT\n";
      }
      if constexpr (support_map_neighbor) {
        std::cout << "supports basic map\n";
      } else {
        std::cerr << "doesn't support basic map, I don't even know how it managed to compile\n";
      }
      if constexpr (support_parallel_map_neighbors) {
        std::cout << "parallel_map_neighbors is supported\n";
      }
      if constexpr (support_map_neighbors_early_exit) {
        std::cout << "map_neighbors_early_exit is supported\n";
      }
      if constexpr (support_parallel_map_neighbors_early_exit) {
        std::cout << "parallel_map_neighbors_early_exit is supported\n";
      }
    } else {
      std::cout << "the graph is directed\n";
      if (support_native_degree_directred) {
        std::cout << "degree is supported\n";
      } else if constexpr (storing_directed_degrees) {
        std::cout << "degrees are stored\n";
      } else {
        std::cerr << "degrees will be calculated as needed\nTHIS IS VERY "
                     "INEFFICIENT\n";
      }
      if constexpr (support_map_in_neighbor && support_map_out_neighbor) {
        std::cout << "supports basic map\n";
      } else {
        std::cerr << "doesn't support basic map, I don't even know how it managed to compile\n";
      }
      if constexpr (support_parallel_map_in_neighbors &&
                    support_parallel_map_out_neighbors) {
        std::cout << "parallel_map_neighbors is supported\n";
      }
      if constexpr (support_map_in_neighbors_early_exit &&
                    support_map_out_neighbors_early_exit) {
        std::cout << "map_neighbors_early_exit is supported\n";
      }
      if constexpr (support_parallel_map_in_neighbors_early_exit &&
                    support_parallel_map_out_neighbors_early_exit) {
        std::cout << "parallel_map_neighbors_early_exit is supported\n";
      }
    }

    if constexpr (support_insert_batch) {
      std::cout << "supports insert_batch\n";
    } else {
      std::cout << "static structure does not support modification\n";
    }
  }

  Graph(Representation &&rep) : g(std::move(rep)) {}

  using weight_type = typename Representation::weight_type;

  size_t N() const {
    if constexpr (n_member) {
      return g.n;
    }
    if constexpr (n_function) {
      return g.N();
    }
  }

  uintE degree(size_t i) const {
    static_assert(symmetric);
    if constexpr (storing_degrees) {
      return stored_degrees[i];
    } else {
      return degree_no_store(i);
    }
  }

  uintE out_degree(size_t i) const {
    if constexpr (support_native_degree_directred && api::degree) {
      return g.out_degree(i);
    } else if constexpr (storing_directed_degrees) {
      return stored_out_degrees[i];
    } else {
      if constexpr (symmetric) {
        return degree(i);
      } else {
        // TODO(wheatman) if we know that this is serial than we don't need a
        // reducer here
        Reducer_sum<size_t> count;
        g.map_out_neighbors(
            i, [&count]([[maybe_unused]] auto... args) { count += 1; });
        return count;
      }
    }
  }

  uintE in_degree(size_t i) const {
    if constexpr (support_native_degree_directred && api::degree) {
      return g.in_degree(i);
    } else if constexpr (storing_directed_degrees) {
      return stored_in_degrees[i];
    } else {
      if constexpr (symmetric) {
        return degree(i);
      } else {
        // TODO(wheatman) if we know that this is serial than we don't need a
        // reducer here
        Reducer_sum<size_t> count;
        g.map_in_neighbors(
            i, [&count]([[maybe_unused]] auto... args) { count += 1; });
        return count;
      }
    }
  }

  void finilize() {
    if constexpr (storing_degrees) {
      stored_degrees.resize(N());
    }
    if constexpr (storing_directed_degrees) {
      stored_in_degrees.resize(N());
      stored_out_degrees.resize(N());
    }
    Reducer_sum<size_t> num_edges;
    parlay::parallel_for(0, N(), [&num_edges, this](auto i) {
      if constexpr (symmetric) {
        uintE deg = degree_no_store(i);
        num_edges += deg;
        if constexpr (storing_degrees) {
          stored_degrees[i] = deg;
        }
      } else {
        uintE in_deg = in_degree_no_store(i);
        uintE out_deg = out_degree_no_store(i);
        num_edges += in_deg + out_deg;
        if constexpr (storing_degrees) {
          stored_in_degrees[i] = in_deg;
          stored_out_degrees[i] = out_deg;
        }
      }
    });
    if constexpr (storing_m) {
      stored_m = num_edges;
    }
  }

  template <class... Args> Graph(Args... args) : g(args...) {
    // static_assert(support_map_neighbor, "must support map\n");
    // static_assert(support_map_in_neighbor || symmetric, "must support map\n");
    // static_assert(support_map_out_neighbor || symmetric, "must support map\n");
    if constexpr (needs_finalize) {
      finilize();
    }
    print_api();
  }

  size_t M() const {
    // static_assert(m_member || m_function, "the underlying graph
    // representation must support getting the number of edges, \n");
    if constexpr (api::M) {
      if constexpr (m_member) {
        return g.m;
      }
      if constexpr (m_function) {
        return g.M();
      }
    }
    if constexpr (storing_m) {
      return stored_m;
    }
    // TODO(wheatman) somehow give a warning in this case since this will be a
    // lot slower
    Reducer_sum<size_t> num_edges;
    parlay::parallel_for(0, N(), [&num_edges, this](auto i) {
      if constexpr (symmetric) {
        num_edges += degree(i);
      } else {
        num_edges += out_degree(i) + in_degree(i);
      }
    });
    return num_edges;
  }

  template <class F> void map_neighbors(size_t id, F f) const {
    static_assert(symmetric);
    g.map_neighbors(id, f);
  }

  template <class F> void map_in_neighbors(size_t id, F f) const {
    if constexpr (support_map_in_neighbor) {
      g.map_in_neighbors(id, f);
    } else {
      map_neighbors(
          id, [f](auto &a, auto &b, auto... args) { return f(b, a, args...); });
    }
  }

  template <class F> void map_out_neighbors(size_t id, F f) const {
    if constexpr (support_map_out_neighbor) {
      g.map_out_neighbors(id, f);
    } else {
      map_neighbors(id, f);
    }
  }

  template <class F> void parallel_map_neighbors(size_t id, F f) const {
    static_assert(symmetric);
    if constexpr (support_parallel_map_neighbors && api::parallel_map) {
      g.parallel_map_neighbors(id, f);
    } else {
      map_neighbors(id, f);
    }
  }

  template <class F> void parallel_map_in_neighbors(size_t id, F f) const {
    if constexpr (support_parallel_map_in_neighbors && api::parallel_map) {
      g.parallel_map_in_neighbors(id, f);
    } else if constexpr (symmetric && api::parallel_map) {
      parallel_map_neighbors(
          id, [f](auto &a, auto &b, auto... args) { return f(b, a, args...); });
    } else {
      map_in_neighbors(id, f);
    }
  }

  template <class F> void parallel_map_out_neighbors(size_t id, F f) const {

    if constexpr (support_parallel_map_out_neighbors && api::parallel_map) {
      g.parallel_map_out_neighbors(id, f);
    } else if constexpr (symmetric && api::parallel_map) {
      parallel_map_neighbors(id, f);
    } else {
      map_out_neighbors(id, f);
    }
  }

  bool parallel_map_desired(size_t id) const {
    // TODO(wheatman) check if it even supports parallel map
    static_assert(symmetric);
    bool desired = api::parallel_map;

    if constexpr (support_native_degree || storing_degrees) {
      if (degree(id) < 1000) {
        desired = false;
      }
    }
    return desired;
  }

  bool parallel_map_in_desired(size_t id) const {
    // TODO(wheatman) check if it even supports parallel map
    bool desired = api::parallel_map;
    if constexpr (support_native_in_degree ||
                  (support_native_degree && symmetric) ||
                  storing_directed_degrees ||
                  (storing_degrees && symmetric )) {
      if (in_degree(id) < 1000) {
        desired = false;
      }
    }
    return desired;
  }

  bool parallel_map_out_desired(size_t id) const {
    // TODO(wheatman) check if it even supports parallel map
    bool desired = api::parallel_map;
    if constexpr (support_native_out_degree ||
                  (support_native_degree && symmetric) ||
                   storing_directed_degrees ||
                  (storing_degrees && symmetric )) {
      if (out_degree(id) < 1000) {
        desired = false;
      }
    }
    return desired;
  }

  template <class F>
  size_t count_neighbors(size_t id, F f, bool parallel = true) const {
    static_assert(symmetric);

    if (parallel_map_desired(id)) {
      Reducer_sum<size_t> count;
      parallel_map_neighbors(id, [&](auto &...args) { count += f(args...); });
      return count;
    } else {
      size_t count = 0;
      map_neighbors(id, [&](auto &...args) { count += f(args...); });
      return count;
    }
  }

  template <class F>
  size_t count_in_neighbors(size_t id, F f, bool parallel = true) const {
    if (parallel_map_in_desired(id)) {
      Reducer_sum<size_t> count;
      parallel_map_in_neighbors(id,
                                [&](auto &...args) { count += f(args...); });
      return count;
    } else {
      size_t count = 0;
      map_in_neighbors(id, [&](auto &...args) { count += f(args...); });
      return count;
    }
  }

  template <class F>
  size_t count_out_neighbors(size_t id, F f, bool parallel = true) const {
    if (parallel_map_out_desired(id) && parallel) {
      Reducer_sum<size_t> count;
      parallel_map_out_neighbors(id,
                                 [&](auto &...args) { count += f(args...); });
      return count;
    } else {
      size_t count = 0;
      map_out_neighbors(id, [&](auto &...args) { count += f(args...); });
      return count;
    }
  }

  template <class F, class Monoid>
  auto map_reduce_neighbors(size_t id, F f, Monoid reduce) const {
    static_assert(symmetric);
    if (parallel_map_desired(id)) {
      Reducer_with_object<Monoid> value(reduce);
      parallel_map_neighbors(id,
                             [&](auto &...args) { value.update(f(args...)); });
      return value.get();
    } else {
      typename Monoid::T value = reduce.identity;
      map_neighbors(
          id, [&](auto &...args) { value = reduce.f(value, f(args...)); });
      return value;
    }
  }

  template <class F, class Monoid>
  auto map_reduce_in_neighbors(size_t id, F f, Monoid reduce) const {
    if (parallel_map_out_desired(id)) {
      Reducer_with_object<Monoid> value(reduce);
      parallel_map_in_neighbors(
          id, [&](auto &...args) { value.update(f(args...)); });
      return value.get();
    } else {
      typename Monoid::T value = reduce.identity;
      map_in_neighbors(
          id, [&](auto &...args) { value = reduce.f(value, f(args...)); });
      return value;
    }
  }

  template <class F, class Monoid>
  auto map_reduce_out_neighbors(size_t id, F f, Monoid reduce) const {
    if (parallel_map_out_desired(id)) {
      Reducer_with_object<Monoid> value(reduce);
      parallel_map_out_neighbors(
          id, [&](auto &...args) { value.update(f(args...)); });
      return value.get();
    } else {
      typename Monoid::T value = reduce.identity;
      map_out_neighbors(
          id, [&](auto &...args) { value = reduce.f(value, f(args...)); });
      return value;
    }
  }

  template <class F> void map_neighbors_early_exit(size_t id, F f) const {
    static_assert(symmetric);

    if constexpr (support_map_neighbors_early_exit && api::map_early_exit) {
      g.map_neighbors_early_exit(id, f);
    } else {
      map_neighbors(id, f);
    }
  }

  template <class F> void map_in_neighbors_early_exit(size_t id, F f) const {

    if constexpr (support_map_in_neighbors_early_exit && api::map_early_exit) {
      g.map_in_neighbors_early_exit(id, f);
    } else {
      if constexpr (symmetric) {
        map_neighbors_early_exit(id, [f](auto &a, auto &b, auto... args) {
          return f(b, a, args...);
        });
      } else {
        map_in_neighbors(id, f);
      }
    }
  }

  template <class F> void map_out_neighbors_early_exit(size_t id, F f) const {

    if constexpr (support_map_out_neighbors_early_exit && api::map_early_exit) {
      g.map_out_neighbors_early_exit(id, f);
    } else {
      if constexpr (symmetric) {
        map_neighbors_early_exit(id, f);
      } else {
        map_out_neighbors(id, f);
      }
    }
  }

  template <class F, typename Block_F = std::nullptr_t>
  void parallel_map_neighbors_early_exit(size_t id, F f,
                                         Block_F block_check = {}) const {
    static_assert(symmetric);

    if constexpr (support_parallel_map_neighbors_early_exit &&
                  api::parallel_map_early_exit) {
      if constexpr (support_parallel_map_neighbors_early_exit1) {
        g.parallel_map_neighbors_early_exit(id, f);
      } else {
        g.parallel_map_neighbors_early_exit(id, f, block_check);
      }
    } else {
      // TODO(wheatman) do we fall back on map_neighbors_early_exit or
      // parallel_map_neighbors
      map_neighbors_early_exit(id, f);
    }
  }

  template <class F, typename Block_F = std::nullptr_t>
  void parallel_map_in_neighbors_early_exit(size_t id, F f,
                                            Block_F block_check = {}) const {

    if constexpr (support_parallel_map_in_neighbors_early_exit &&
                  api::parallel_map_early_exit) {
      if constexpr (support_parallel_map_in_neighbors_early_exit1) {
        g.parallel_map_in_neighbors_early_exit(id, f);
      } else {
        g.parallel_map_in_neighbors_early_exit(id, f, block_check);
      }

    } else {
      if constexpr (symmetric) {
        parallel_map_neighbors_early_exit(
            id,
            [f](auto &a, auto &b, auto... args) { return f(b, a, args...); },
            block_check);
      } else {
        // TODO(wheatman) do we fall back on map_neighbors_early_exit or
        // parallel_map_neighbors
        map_in_neighbors_early_exit(id, f);
      }
    }
  }

  template <class F, typename Block_F = std::nullptr_t>
  void parallel_map_out_neighbors_early_exit(size_t id, F f,
                                             Block_F block_check = {}) const {

    if constexpr (support_parallel_map_out_neighbors_early_exit &&
                  api::parallel_map_early_exit) {
      if constexpr (support_parallel_map_out_neighbors_early_exit1) {
        g.parallel_map_out_neighbors_early_exit(id, f);
      } else {
        g.parallel_map_out_neighbors_early_exit(id, f, block_check);
      }

    } else {
      if constexpr (symmetric) {
        parallel_map_neighbors_early_exit(id, f, block_check);
      } else {
        // TODO(wheatman) do we fall back on map_neighbors_early_exit or
        // parallel_map_neighbors
        map_out_neighbors_early_exit(id, f);
      }
    }
  }
#if false
  //TODO(wheatman) add the optimization for these guys in
  template <class F>
  void map_neighbor_range(size_t start_id, size_t end_id, F f) const {
    static_assert(symmetric);
    constexpr bool support_native =
        requires(const Representation &g, size_t start_i, size_t end_i, F f) {
      g.map_neighbor_range(start_i, end_i, f);
    };
    if constexpr (support_native && api::map_range) {
      g.map_neighbor_range(start_id, end_id, f);
    } else {
      for (size_t i = start_id; i < end_id; i++) {
        map_neighbors(i, f);
      }
    }
  }

  template <class F>
  void map_in_neighbor_range(size_t start_id, size_t end_id, F f) const {
    constexpr bool support_native =
        requires(const Representation &g, size_t start_i, size_t end_i, F f) {
      g.map_in_neighbor_range(start_i, end_i, f);
    };
    if constexpr (support_native && api::map_range) {
      g.map_in_neighbor_range(start_id, end_id, f);
    } else {
      if constexpr (symmetric) {
        map_neighbor_range(
            start_id, end_id,
            [f](auto &a, auto &b, auto... args) { return f(b, a, args...); });
      } else {
        for (size_t i = start_id; i < end_id; i++) {
          map_in_neighbors(i, f);
        }
      }
    }
  }

  template <class F>
  void map_out_neighbor_range(size_t start_id, size_t end_id, F f) const {
    constexpr bool support_native =
        requires(const Representation &g, size_t start_i, size_t end_i, F f) {
      g.map_out_neighbor_range(start_i, end_i, f);
    };
    if constexpr (support_native && api::map_range) {
      g.map_out_neighbor_range(start_id, end_id, f);
    } else {
      if constexpr (symmetric) {
        map_neighbor_range(start_id, end_id, f);
      } else {
        for (size_t i = start_id; i < end_id; i++) {
          map_out_neighbors(i, f);
        }
      }
    }
  }
#endif

sequence<std::tuple<uintE, uintE, weight_type>> edges() const {
    static_assert(symmetric);
    using g_edge = std::tuple<uintE, uintE, weight_type>;
    auto degs = sequence<size_t>::from_function(
        N(), [&](size_t i) { return out_degree(i); });
    size_t sum_degs = parlay::scan_inplace(make_slice(degs));
    assert(sum_degs == M());
    auto edges = sequence<g_edge>(sum_degs);
    parallel_for(0, N(),
                 [&](size_t i) {
                   size_t k = degs[i];
                   auto map_f = [&](const uintE& u, const uintE& v,
                                    const weight_type& wgh) {
                     edges[k++] = std::make_tuple(u, v, wgh);
                   };
                   map_neighbors(i, map_f);
                 },
                 1);
    return edges;
  }

  void print_adj() const {
    if constexpr (needs_finalize) {
      finilize();
    }
    std::cout << "AdjacencyGraph\n";
    std::cout << N() << "\n";
    std::cout << M() << "\n";
    size_t running_count = 0;;
    for (uint64_t i = 0; i < N(); i++) {
      std::cout << running_count << "\n";
      running_count += degree(i);
    }
    for (uint64_t i = 0; i < N(); i++) {
      map_neighbors(i, [&](auto src, auto dest, auto val) {
        std::cout << dest << "\n";
      });
    }
  }
  void write_adj(std::string fname) {
    if constexpr (needs_finalize) {
      finilize();
    }
    std::ofstream myfile;
    myfile.open(fname);
        
    myfile << "AdjacencyGraph\n";
    myfile << N() << "\n";
    myfile << M() << "\n";
    size_t running_count = 0;;
    for (uint64_t i = 0; i < N(); i++) {
      myfile << running_count << "\n";
      running_count += degree(i);
    }
    for (uint64_t i = 0; i < N(); i++) {
      map_neighbors(i, [&](auto src, auto dest, auto val) {
        myfile << dest << "\n";
      });
    }
    myfile.close();
  }

  void insert_sorted_batch(std::tuple<gbbs::uintE, gbbs::uintE> *es, size_t n) {
    g.insert_sorted_batch(es, n);
  }

  void remove_sorted_batch(std::tuple<gbbs::uintE, gbbs::uintE> *es, size_t n) {
    g.remove_sorted_batch(es, n);
  }

  void insert_sorted_grouped_batch(parlay::sequence<std::pair<gbbs::uintE, parlay::sequence<gbbs::uintE>>> &els) {
    g.insert_sorted_grouped_batch(els);

  }

  void remove_sorted_grouped_batch(parlay::sequence<std::pair<gbbs::uintE, parlay::sequence<gbbs::uintE>>> & els) {
    g.remove_sorted_grouped_batch(els);
  }

  void insert_batch(std::tuple<gbbs::uintE, gbbs::uintE> *es, size_t n) {
    if constexpr (support_insert_grouped_batch) {
      auto groups = semisort::group_by(parlay::slice(es, es+n), 
        [](auto elem) {return std::get<0>(elem);}, 
        [](auto elem) {return std::get<1>(elem);});
      parlay::parallel_for(0, groups.size(), [&](size_t i) {
        parlay::integer_sort_inplace(groups[i].second);
        auto seq = parlay::unique(groups[i].second);
        groups[i].second = seq;
      });
      insert_sorted_grouped_batch(groups);
    } else {
      auto elements = parlay::unique(parlay::sort(parlay::slice(es, es+n)));
      insert_sorted_batch(elements.data(), elements.size());

    }
  }
  void remove_batch(std::tuple<gbbs::uintE, gbbs::uintE> *es, size_t n) {
    if constexpr (support_insert_grouped_batch) {
      auto groups = semisort::group_by(parlay::slice(es, es+n), 
        [](auto elem) {return std::get<0>(elem);}, 
        [](auto elem) {return std::get<1>(elem);});
      parlay::parallel_for(0, groups.size(), [&](size_t i) {
        parlay::integer_sort_inplace(groups[i].second);
        auto seq = parlay::unique(groups[i].second);
        groups[i].second = seq;
      });
      remove_sorted_grouped_batch(groups);
    } else {
      auto elements = parlay::unique(parlay::sort(parlay::slice(es, es+n)));
      remove_sorted_batch(elements.data(), elements.size());

    }
  }
};

// TODO(wheatman) make sure that these are only being used for intermediate
// objects, that it is ok that we control the API for
// maybe add the ability to specify which represenation to use?
template <class graph, class EdgeSeq, class GetU, class GetV, class GetW>
static inline graph
sym_graph_from_edges(EdgeSeq &A, size_t n, GetU &&get_u, GetV &&get_v,
                     GetW &&get_w, bool is_sorted = false) {
  return graph_implementations::template sym_graph_from_edges<graph>(A, n, get_u, get_v, get_w,
                                                     is_sorted);
}

template <class graph>
static inline graph
sym_graph_from_edges(sequence<std::tuple<uintE, uintE, typename graph::weight_type>> &A, size_t n,
                     bool is_sorted = false) {
  return graph_implementations::template sym_graph_from_edges<graph>(A, n, is_sorted);
}

}  // namespace gbbs

#ifdef GRAPH_API_DEGREE
constexpr bool graph_api_use_degree = GRAPH_API_DEGREE;
#else
constexpr bool graph_api_use_degree = true;
#endif
#ifdef GRAPH_API_M
constexpr bool graph_api_M = GRAPH_API_M;
#else
constexpr bool graph_api_M = true;
#endif
#ifdef GRAPH_API_PARALLEL_MAP
constexpr bool graph_api_parallel_map = GRAPH_API_PARALLEL_MAP;
#else
constexpr bool graph_api_parallel_map = true;
#endif
#ifdef GRAPH_API_MAP_EARLY_EXIT
constexpr bool graph_api_map_early_exit = GRAPH_API_MAP_EARLY_EXIT;
#else
constexpr bool graph_api_map_early_exit = true;
#endif
#ifdef GRAPH_API_PARALLEL_MAP_EARLY_EXIT
constexpr bool graph_api_parallel_map_early_exit =
    GRAPH_API_PARALLEL_MAP_EARLY_EXIT;
#else
constexpr bool graph_api_parallel_map_early_exit = true;
#endif
#ifdef GRAPH_API_MAP_RANGE
constexpr bool graph_api_map_range = GRAPH_API_MAP_RANGE;
#else
constexpr bool graph_api_map_range = true;
#endif
#ifdef GRAPH_API_STORE_M
constexpr bool graph_api_store_m = GRAPH_API_STORE_M;
#else
constexpr bool graph_api_store_m = true;
#endif
#ifdef GRAPH_API_STORE_DEGREES
constexpr bool graph_api_store_degrees = GRAPH_API_STORE_DEGREES;
#else
constexpr bool graph_api_store_degrees = true;
#endif

#ifdef UNWEIGHTED_SYM_GRAPH_IMPL
using unweighted_sym_graph = gbbs::Graph<UNWEIGHTED_SYM_GRAPH_IMPL, true, gbbs::Graph_API_use<graph_api_use_degree, graph_api_M,
                        graph_api_parallel_map, graph_api_map_early_exit,
                        graph_api_parallel_map_early_exit, graph_api_map_range, graph_api_store_m, graph_api_store_degrees>>;
#else
using unweighted_sym_graph = gbbs::Graph<gbbs::graph_implementations::symmetric_graph<
      gbbs::symmetric_vertex, gbbs::empty>, true, gbbs::Graph_API_use<graph_api_use_degree, graph_api_M,
                        graph_api_parallel_map, graph_api_map_early_exit,
                        graph_api_parallel_map_early_exit, graph_api_map_range, graph_api_store_m, graph_api_store_degrees>>;
#endif


using unweighted_asym_graph = gbbs::Graph<gbbs::graph_implementations::asymmetric_graph<gbbs::asymmetric_vertex, gbbs::empty>, false>;

template <class Weight>
using weighted_sym_graph = gbbs::Graph<gbbs::graph_implementations::symmetric_graph<gbbs::symmetric_vertex, Weight>, true, gbbs::full_api>;
template <class Weight>
using weighted_asym_graph = gbbs::Graph<gbbs::graph_implementations::asymmetric_graph<gbbs::asymmetric_vertex, Weight>, false>;
