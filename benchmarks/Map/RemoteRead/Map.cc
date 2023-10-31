#include "Map.h"

namespace gbbs {

template <class Graph> double Map_runner(Graph &G, commandLine P) {
  std::cout << "### Application: Test Map" << std::endl;
  size_t read_size = static_cast<uintE>(P.getOptionLongValue("-rs", 0));
  // read size of 0 means 1 bit, else in a bumber of bytes
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.N() << std::endl;
  std::cout << "### m: " << G.M() << std::endl;
  std::cout << "### ------------------------------------" << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  size_t n = G.N();
  if (read_size == 0) {
    parlay::random r = random();
    std::vector<bool> remote_data(n);
    remote_data[r.rand() % n] = true;
    timer t;
    t.start();
    auto sum = Test_Map_Remote(G, remote_data);
    double tt = t.stop();
    std::cout << "edge sum was " << sum << std::endl;

    std::cout << "### Running Time: " << tt << std::endl;
    return tt;
  }
  if (read_size == 1) {
    parlay::random r = random();
    std::vector<char> remote_data(n);
    remote_data[r.rand() % n] = 1;
    timer t;
    t.start();
    auto sum = Test_Map_Remote(G, remote_data);
    double tt = t.stop();
    std::cout << "edge sum was " << sum << std::endl;

    std::cout << "### Running Time: " << tt << std::endl;
    return tt;
  }
  if (read_size == 2) {
    parlay::random r = random();
    std::vector<short> remote_data(n);
    remote_data[r.rand() % n] = 1;
    timer t;
    t.start();
    auto sum = Test_Map_Remote(G, remote_data);
    double tt = t.stop();
    std::cout << "edge sum was " << sum << std::endl;

    std::cout << "### Running Time: " << tt << std::endl;
    return tt;
  }
  if (read_size == 4) {
    parlay::random r = random();
    std::vector<int> remote_data(n);
    remote_data[r.rand() % n] = 1;
    timer t;
    t.start();
    auto sum = Test_Map_Remote(G, remote_data);
    double tt = t.stop();
    std::cout << "edge sum was " << sum << std::endl;

    std::cout << "### Running Time: " << tt << std::endl;
    return tt;
  }
  if (read_size == 8) {
    parlay::random r = random();
    std::vector<long> remote_data(n);
    remote_data[r.rand() % n] = 1;
    timer t;
    t.start();
    auto sum = Test_Map_Remote(G, remote_data);
    double tt = t.stop();
    std::cout << "edge sum was " << sum << std::endl;

    std::cout << "### Running Time: " << tt << std::endl;
    return tt;
  }
  if (read_size == 16) {
    parlay::random r = random();
    std::vector<wide_int<size_t, 16>> remote_data(n);
    remote_data[r.rand() % n] = 1;
    timer t;
    t.start();
    auto sum = Test_Map_Remote(G, remote_data);
    double tt = t.stop();
    std::cout << "edge sum was " << sum << std::endl;

    std::cout << "### Running Time: " << tt << std::endl;
    return tt;
  }
  if (read_size == 32) {
    parlay::random r = random();
    std::vector<wide_int<size_t, 16>> remote_data(n);
    remote_data[r.rand() % n] = 1;
    timer t;
    t.start();
    auto sum = Test_Map_Remote(G, remote_data);
    double tt = t.stop();
    std::cout << "edge sum was " << sum << std::endl;

    std::cout << "### Running Time: " << tt << std::endl;
    return tt;
  }
  if (read_size == 64) {
    parlay::random r = random();
    std::vector<wide_int<size_t, 16>> remote_data(n);
    remote_data[r.rand() % n] = 1;
    timer t;
    t.start();
    auto sum = Test_Map_Remote(G, remote_data);
    double tt = t.stop();
    std::cout << "edge sum was " << sum << std::endl;

    std::cout << "### Running Time: " << tt << std::endl;
    return tt;
  }
  return 0;
}
} // namespace gbbs

generate_symmetric_main(gbbs::Map_runner, false);