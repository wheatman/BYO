#include "gbbs/graph_io.h"

#include <vector>

#include "gbbs/unit_tests/graph_test_utils.h"
#include "gtest/gtest.h"

namespace gbbs {

namespace gi = gbbs_io;
namespace gt = graph_test;
using NoWeight = gbbs::empty;

TEST(EdgeListToAsymmetricGraph, NoEdges) {
  const std::vector<gi::Edge<NoWeight>> kEdges{};
  const auto graph{gi::edge_list_to_asymmetric_graph(kEdges)};
  EXPECT_EQ(graph.N(), 0);
  EXPECT_EQ(graph.M(), 0);
}

TEST(EdgeListToAsymmetricGraph, DuplicateEdges) {
  // Check that `EdgeListToAsymmetricGraph()` works even when there are
  // duplicate edges in the input.
  //
  // Graph diagram:
  // 0 --> 1
  const std::vector<gi::Edge<NoWeight>> kEdges{
      {0, 1}, {0, 1},
  };
  auto graph{gi::edge_list_to_asymmetric_graph(kEdges)};
  EXPECT_EQ(graph.N(), 2);
  EXPECT_EQ(graph.M(), 1);

  {
    const std::vector<uintE> kExpectedOutNeighbors{1};
    gt::CheckUnweightedOutNeighbors(graph, 0, kExpectedOutNeighbors);
    EXPECT_EQ(graph.in_degree(0), 0);
  }
  {
    const std::vector<uintE> kExpectedInNeighbors{0};
    EXPECT_EQ(graph.out_degree(1), 0);
    gt::CheckUnweightedInNeighbors(graph, 1, kExpectedInNeighbors);
  }
}

TEST(EdgeListToAsymmetricGraph, SkipFirstVertex) {
  // Check that `EdgeListToAsymmetricGraph()` works even when the first vertex
  // has degree 0.
  //
  // Graph diagram:
  // 1 --> 2
  const std::vector<gi::Edge<NoWeight>> kEdges{
      {1, 2},
  };
  auto graph{gi::edge_list_to_asymmetric_graph(kEdges)};
  EXPECT_EQ(graph.N(), 3);
  EXPECT_EQ(graph.M(), kEdges.size());

  {
    EXPECT_EQ(graph.out_degree(0), 0);
    EXPECT_EQ(graph.in_degree(0), 0);
  }
  {
    const std::vector<uintE> kExpectedOutNeighbors{2};
    gt::CheckUnweightedOutNeighbors(graph, 1, kExpectedOutNeighbors);
    EXPECT_EQ(graph.in_degree(1), 0);
  }
  {
    const std::vector<uintE> kExpectedInNeighbors{1};
    EXPECT_EQ(graph.out_degree(2), 0);
    gt::CheckUnweightedInNeighbors(graph, 2, kExpectedInNeighbors);
  }
}

TEST(EdgeListToAsymmetricGraph, OutOfOrderEdges) {
  // Check that `EdgeListToAsymmetricGraph()` works even when the input edge
  // list is in a scrambled order.
  //
  // Graph diagram:
  // 0 --> 1     3 --> 6
  // | ^   ^
  // |  \  |     4
  // v   \ |
  // 2 <-> 5

  const std::vector<gi::Edge<NoWeight>> kEdges{
      {3, 6}, {0, 2}, {5, 0}, {5, 1}, {2, 5}, {0, 1}, {5, 2},
  };
  auto graph{gi::edge_list_to_asymmetric_graph(kEdges)};
  EXPECT_EQ(graph.N(), 7);
  EXPECT_EQ(graph.M(), kEdges.size());

  {
    const std::vector<uintE> kExpectedOutNeighbors{1, 2};
    const std::vector<uintE> kExpectedInNeighbors{5};
    gt::CheckUnweightedOutNeighbors(graph, 0, kExpectedOutNeighbors);
    gt::CheckUnweightedInNeighbors(graph, 0, kExpectedInNeighbors);
  }
  {
    const std::vector<uintE> kExpectedOutNeighbors{};
    const std::vector<uintE> kExpectedInNeighbors{0, 5};
    gt::CheckUnweightedOutNeighbors(graph, 1, kExpectedOutNeighbors);
    gt::CheckUnweightedInNeighbors(graph, 1, kExpectedInNeighbors);
  }
  {
    const std::vector<uintE> kExpectedOutNeighbors{5};
    const std::vector<uintE> kExpectedInNeighbors{0, 5};
    gt::CheckUnweightedOutNeighbors(graph, 2, kExpectedOutNeighbors);
    gt::CheckUnweightedInNeighbors(graph, 2, kExpectedInNeighbors);
  }
  {
    const std::vector<uintE> kExpectedOutNeighbors{6};
    const std::vector<uintE> kExpectedInNeighbors{};
    gt::CheckUnweightedOutNeighbors(graph, 3, kExpectedOutNeighbors);
    gt::CheckUnweightedInNeighbors(graph, 3, kExpectedInNeighbors);
  }
  {
    const std::vector<uintE> kExpectedOutNeighbors{};
    const std::vector<uintE> kExpectedInNeighbors{};
    gt::CheckUnweightedOutNeighbors(graph, 4, kExpectedOutNeighbors);
    gt::CheckUnweightedInNeighbors(graph, 4, kExpectedInNeighbors);
  }
  {
    const std::vector<uintE> kExpectedOutNeighbors{0, 1, 2};
    const std::vector<uintE> kExpectedInNeighbors{2};
    gt::CheckUnweightedOutNeighbors(graph, 5, kExpectedOutNeighbors);
    gt::CheckUnweightedInNeighbors(graph, 5, kExpectedInNeighbors);
  }
  {
    const std::vector<uintE> kExpectedOutNeighbors{};
    const std::vector<uintE> kExpectedInNeighbors{3};
    gt::CheckUnweightedOutNeighbors(graph, 6, kExpectedOutNeighbors);
    gt::CheckUnweightedInNeighbors(graph, 6, kExpectedInNeighbors);
  }
}

TEST(EdgeListToSymmetricGraph, NoEdges) {
  const std::vector<gi::Edge<NoWeight>> kEdges{};
  const auto graph{gi::edge_list_to_symmetric_graph(kEdges)};
  EXPECT_EQ(graph.N(), 0);
  EXPECT_EQ(graph.M(), 0);
}

TEST(EdgeListToSymmetricGraph, DuplicateEdges) {
  // Check that `EdgeListToSymmetricGraph()` works even when there are
  // duplicate edges in the input.
  //
  // Graph diagram:
  // 0 --- 1
  const std::vector<gi::Edge<NoWeight>> kEdges{
      {0, 1}, {1, 0}, {0, 1},
  };
  auto graph{gi::edge_list_to_symmetric_graph(kEdges)};
  EXPECT_EQ(graph.N(), 2);
  EXPECT_EQ(graph.M(), 2);

  {
    const std::vector<uintE> kExpectedOutNeighbors{1};
    gt::CheckUnweightedOutNeighbors(graph, 0, kExpectedOutNeighbors);
  }
  {
    const std::vector<uintE> kExpectedOutNeighbors{0};
    gt::CheckUnweightedOutNeighbors(graph, 1, kExpectedOutNeighbors);
  }
}

TEST(EdgeListToSymmetricGraph, SkipFirstVertex) {
  // Check that `EdgeListToSymmetricGraph()` works even when the first vertex
  // has degree 0.
  //
  // Graph diagram:
  // 1 --- 2
  const std::vector<gi::Edge<NoWeight>> kEdges{
      {1, 2},
  };
  auto graph{gi::edge_list_to_symmetric_graph(kEdges)};
  EXPECT_EQ(graph.N(), 3);
  EXPECT_EQ(graph.M(), 2);

  {
    EXPECT_EQ(graph.out_degree(0), 0);
  }
  {
    const std::vector<uintE> kExpectedOutNeighbors{2};
    gt::CheckUnweightedOutNeighbors(graph, 1, kExpectedOutNeighbors);
  }
  {
    const std::vector<uintE> kExpectedOutNeighbors{1};
    gt::CheckUnweightedOutNeighbors(graph, 2, kExpectedOutNeighbors);
  }
}

TEST(EdgeListToSymmetricGraph, OutOfOrderEdges) {
  // Check that `EdgeListToSymmetricGraph()` works even when the input edge
  // list is in a scrambled order.
  //
  // Graph diagram:
  // 0 --- 1     3 --- 6
  // | \   |
  // |  \  |     4
  // |   \ |
  // 2     5

  const std::vector<gi::Edge<NoWeight>> kEdges{
      {1, 5}, {0, 5}, {6, 3}, {2, 0}, {1, 0},
  };
  auto graph{gi::edge_list_to_symmetric_graph(kEdges)};
  EXPECT_EQ(graph.N(), 7);
  EXPECT_EQ(graph.M(), 10);

  {
    const std::vector<uintE> kExpectedOutNeighbors{1, 2, 5};
    gt::CheckUnweightedOutNeighbors(graph, 0, kExpectedOutNeighbors);
  }
  {
    const std::vector<uintE> kExpectedOutNeighbors{0, 5};
    gt::CheckUnweightedOutNeighbors(graph, 1, kExpectedOutNeighbors);
  }
  {
    const std::vector<uintE> kExpectedOutNeighbors{0};
    gt::CheckUnweightedOutNeighbors(graph, 2, kExpectedOutNeighbors);
  }
  {
    const std::vector<uintE> kExpectedOutNeighbors{6};
    gt::CheckUnweightedOutNeighbors(graph, 3, kExpectedOutNeighbors);
  }
  {
    const std::vector<uintE> kExpectedOutNeighbors{};
    gt::CheckUnweightedOutNeighbors(graph, 4, kExpectedOutNeighbors);
  }
  {
    const std::vector<uintE> kExpectedOutNeighbors{0, 1};
    gt::CheckUnweightedOutNeighbors(graph, 5, kExpectedOutNeighbors);
  }
  {
    const std::vector<uintE> kExpectedOutNeighbors{3};
    gt::CheckUnweightedOutNeighbors(graph, 6, kExpectedOutNeighbors);
  }
}

}  // namespace gbbs
