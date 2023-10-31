// Utility functions for working with graphs that are useful for writing unit
// tests.
#pragma once

#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include "gbbs/bridge.h"
#include "gbbs/graph.h"
#include "gbbs/helpers/undirected_edge.h"
#include "gbbs/macros.h"
#include "gbbs/vertex.h"
#include "gtest/gtest.h"

namespace gbbs {
namespace graph_test {

namespace internal {  // Internal declarations

template <typename Weight>
void CheckUnweightedNeighbors(size_t actual_degree,
                              const std::tuple<uintE, Weight>* actual_neighbors,
                              const std::vector<uintE>& expected_neighbors);

template <typename Weight>
void CheckWeightedNeighbors(
    size_t actual_degree, const std::tuple<uintE, Weight>* actual_neighbors,
    const std::vector<std::tuple<uintE, Weight>>& expected_neighbors);

}  // namespace internal

// Make an undirected, unweighted graph from a list of edges.
template <class graph>
graph MakeUnweightedSymmetricGraph(const uintE num_vertices,
                             const std::unordered_set<UndirectedEdge> &edges) {
  using Edge = std::tuple<uintE, uintE, gbbs::empty>;
  constexpr gbbs::empty weight{};
  sequence<Edge> edge_sequence(edges.size() * 2);
  auto edges_it{edges.cbegin()};
  for (size_t i = 0; i < edges.size(); i++) {
    edge_sequence[2 * i] = std::make_tuple(
        edges_it->endpoints().first, edges_it->endpoints().second, weight);
    edge_sequence[2 * i + 1] = std::make_tuple(
        edges_it->endpoints().second, edges_it->endpoints().first, weight);
    ++edges_it;
  }
  // TODO(tomtseng): Some graph operations assume that the neighbor lists are
  // sorted, so we sort here. But maybe this sorting belongs in
  // `sym_graph_from_edges` or as an option to `sym_graph_from_edges`.
  // See https://github.com/ldhulipala/gbbs/pull/21.
  parlay::sample_sort_inplace(
      make_slice(edge_sequence), [](const Edge &left, const Edge &right) {
        return std::tie(std::get<0>(left), std::get<1>(left)) <
               std::tie(std::get<0>(right), std::get<1>(right));
      });
  constexpr bool kEdgesAreSorted{true};
  return sym_graph_from_edges<graph>(edge_sequence, num_vertices, kEdgesAreSorted);
}

// Check that vertex has `expected_neighbors` as its out-neighbors. Does not
// check edge weights. Ordering matters.
template <class Graph>
void CheckUnweightedOutNeighbors(Graph& g, size_t i,
                                 const std::vector<uintE>& expected_neighbors) {
  EXPECT_EQ(g.out_degree(i), expected_neighbors.size());
  std::vector<uintE> actual_neighbors;
  g.map_out_neighbors(i, [&](auto src, auto dest, auto weight) {actual_neighbors.emplace_back(dest);});
  EXPECT_EQ(actual_neighbors, expected_neighbors);
}

// Check that vertex has `expected_neighbors` as its in-neighbors. Does not
// check edge weights. Ordering matters.
template <class Graph>
void CheckUnweightedInNeighbors(Graph& g, size_t i,
                                const std::vector<uintE>& expected_neighbors) {
  EXPECT_EQ(g.in_degree(i), expected_neighbors.size());
  std::vector<uintE> actual_neighbors;
  g.map_in_neighbors(i, [&](auto src, auto dest, auto weight) {actual_neighbors.emplace_back(src);});
  EXPECT_EQ(actual_neighbors, expected_neighbors);
}

// Check that vertex has `expected_neighbors` as its out-neighbors. Does not
// check edge weights. Ordering matters.
template <class Vertex>
void CheckWeightedOutNeighbors(
    Vertex& vertex,
    const std::vector<std::tuple<uintE, typename Vertex::weight_type>>&
        expected_neighbors) {
  internal::CheckWeightedNeighbors(vertex.out_degree(),
                                   vertex.out_neighbors().neighbors,
                                   expected_neighbors);
}

// Check that vertex has `expected_neighbors` as its in-neighbors. Does not
// check edge weights. Ordering matters.
template <class Vertex>
void CheckWeightedInNeighbors(
    Vertex& vertex,
    const std::vector<std::tuple<uintE, typename Vertex::weight_type>>&
        expected_neighbors) {
  internal::CheckWeightedNeighbors(
      vertex.in_degree(), vertex.in_neighbors().neighbors, expected_neighbors);
}

namespace internal {  // Internal definitions

// Check that the first `actual_degree` entries of `actual_neighbors` equal
// the entire vector `expected_neighbors`.
template <typename Weight>
void CheckWeightedNeighbors(
    const size_t actual_degree,
    const std::tuple<uintE, Weight>* const actual_neighbors,
    const std::vector<std::tuple<uintE, Weight>>& expected_neighbors) {
  EXPECT_EQ(actual_degree, expected_neighbors.size());

  std::vector<std::tuple<uintE, Weight>> neighbors_vector;
  neighbors_vector.reserve(actual_degree);
  for (size_t i = 0; i < actual_degree; i++) {
    neighbors_vector.emplace_back(std::get<0>(actual_neighbors[i]),
                                  std::get<1>(actual_neighbors[i]));
  }

  EXPECT_EQ(neighbors_vector, expected_neighbors);
}

}  // namespace internal

}  // namespace graph_test
}  // namespace gbbs
