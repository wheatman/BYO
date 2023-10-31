#include "Map.h"

namespace gbbs {

template <class Graph> double Map_runner(Graph &G, commandLine P) {
  std::cout << "### Application: Test Map" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.N() << std::endl;
  std::cout << "### m: " << G.M() << std::endl;
  std::cout << "### ------------------------------------" << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  timer t;
  t.start();
  auto sum = Test_Map(G);
  double tt = t.stop();
  std::cout << "edge sum was " << sum << std::endl;

  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}
} // namespace gbbs

generate_symmetric_main(gbbs::Map_runner, false);