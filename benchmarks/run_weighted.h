#pragma once

#include <map>

#include "GeneralWeightSSSP/BellmanFord/BellmanFord.h"
#include "IntegralWeightSSSP/JulienneDBS17/wBFS.h"
#include "PositiveWeightSSSP/DeltaStepping/DeltaStepping.h"
#include "SSWidestPath/JulienneDBS17/SSWidestPath.h"

namespace gbbs {

template <class Graph>
double WBFS_runner(const Graph &G, uintE src, size_t rounds, size_t num_buckets,
                   bool no_blocked, bool largemem, bool dump) {
  std::cout << "### Application: wBFS (Weighted Breadth-First Search)"
            << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.N() << std::endl;
  std::cout << "### m: " << G.M() << std::endl;
  std::cout << "### Params: -src = " << src
            << " -nb (num_buckets) = " << num_buckets << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  if (num_buckets != (((uintE)1) << parlay::log2_up(num_buckets))) {
    std::cout << "Please specify a number of buckets that is a power of two"
              << "\n";
    exit(-1);
  }

  double total_time = 0.0;
  for (size_t r = 0; r <= rounds; r++) {
    timer t;
    t.start();
    auto dists = wBFS(G, src, num_buckets, largemem, no_blocked);
    double tt = t.stop();
    std::cout << "### Running Time: " << tt << std::endl;
    if (r == 0) {
      if (dump) {
        // useful for debugging
        std::ofstream myfile;
        myfile.open("wbfs.out");
        for (unsigned int i = 0; i < G.N(); i++) {
          myfile << dists[i] << std::endl;
        }
        myfile.close();
      }
    } else {
      total_time += tt;
    }
  }
  auto time_per_iter = total_time / rounds;
  std::cout << "# time per iter: " << time_per_iter << "\n";
  return time_per_iter;
}

template <class Graph>
double BellmanFord_runner(const Graph &G, uintE src, size_t rounds, bool dump) {
  std::cout << "### Application: BellmanFord" << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.N() << std::endl;
  std::cout << "### m: " << G.M() << std::endl;
  std::cout << "### Params: -src = " << src << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  double total_time = 0.0;
  for (size_t r = 0; r <= rounds; r++) {
    timer t;
    t.start();
    auto dists = BellmanFord(G, src);
    double tt = t.stop();
    std::cout << "### Running Time: " << tt << std::endl;
    if (r == 0) {
      if (dump) {
        // useful for debugging
        std::ofstream myfile;
        myfile.open("bf.out");
        for (unsigned int i = 0; i < G.N(); i++) {
          myfile << dists[i] << std::endl;
        }
        myfile.close();
      }
    } else {
      total_time += tt;
    }
  }
  auto time_per_iter = total_time / rounds;
  std::cout << "# time per iter: " << time_per_iter << "\n";
  return time_per_iter;
}

template <class Graph>
double DeltaStepping_runner(const Graph &G, uintE src, size_t rounds,
                            double delta, size_t num_buckets, bool dump) {
  std::cout << "### Application: wBFS (Weighted Breadth-First Search)"
            << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.N() << std::endl;
  std::cout << "### m: " << G.M() << std::endl;
  std::cout << "### Params: -src = " << src
            << " -nb (num_buckets) = " << num_buckets << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  if (num_buckets != (((uintE)1) << parlay::log2_up(num_buckets))) {
    std::cout << "Please specify a number of buckets that is a power of two"
              << "\n";
    exit(-1);
  }

  double total_time = 0.0;
  for (size_t r = 0; r <= rounds; r++) {
    timer t;
    t.start();
    auto dists = DeltaStepping(G, src, delta, num_buckets);
    double tt = t.stop();
    std::cout << "### Running Time: " << tt << std::endl;
    if (r == 0) {
      if (dump) {
        // useful for debugging
        std::ofstream myfile;
        myfile.open("ds.out");
        for (unsigned int i = 0; i < G.N(); i++) {
          myfile << dists[i] << std::endl;
        }
        myfile.close();
      }
    } else {
      total_time += tt;
    }
  }
  auto time_per_iter = total_time / rounds;
  std::cout << "# time per iter: " << time_per_iter << "\n";
  return time_per_iter;
}

template <class Graph>
double SSWidestPath_runner(const Graph &G, uintE src, size_t rounds,
                           size_t num_buckets, bool no_blocked, bool largemem,
                           bool dump) {
  std::cout << "### Application: wBFS (Weighted Breadth-First Search)"
            << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.N() << std::endl;
  std::cout << "### m: " << G.M() << std::endl;
  std::cout << "### Params: -src = " << src
            << " -nb (num_buckets) = " << num_buckets << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  if (num_buckets != (((uintE)1) << parlay::log2_up(num_buckets))) {
    std::cout << "Please specify a number of buckets that is a power of two"
              << "\n";
    exit(-1);
  }

  double total_time = 0.0;
  for (size_t r = 0; r <= rounds; r++) {
    timer t;
    t.start();
    auto widths = SSWidestPath(G, src, num_buckets, largemem, no_blocked);

    double tt = t.stop();
    std::cout << "### Running Time: " << tt << std::endl;
    if (r == 0) {
      if (dump) {
        // useful for debugging
        std::ofstream myfile;
        myfile.open("wp.out");
        for (unsigned int i = 0; i < G.N(); i++) {
          myfile << widths[i] << std::endl;
        }
        myfile.close();
      }
    } else {
      total_time += tt;
    }
  }
  auto time_per_iter = total_time / rounds;
  std::cout << "# time per iter: " << time_per_iter << "\n";
  return time_per_iter;
}

class run_all_options {
public:
  uintE src = 0;
  size_t rounds = 3;
  bool dump = false;
  bool verify = false;
  size_t wbfs_num_buckets = 32;
  bool wbfs_no_blocked = false;
  bool wbfs_largemem = false;
  double ds_delta = 1.0;
  size_t ds_num_buckets = 32;
  size_t wp_num_buckets = 32;
  bool wp_no_blocked = false;
  bool wp_largemem = false;
};

template <class Graph>
void run_all(const Graph &G, const run_all_options &options) {
  // static constexpr bool symmetric = Graph::symmetric;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.N() << std::endl;
  std::cout << "### m: " << G.M() << std::endl;
  if (options.dump) {
    std::cout << "writing output arrays to files\n";
  }
  std::map<std::string, double> time_map;

  time_map["WBFS"] =
      WBFS_runner(G, options.src, options.rounds, options.wbfs_num_buckets,
                  options.wbfs_no_blocked, options.wbfs_largemem, options.dump);

  time_map["BellmanFord"] =
      BellmanFord_runner(G, options.src, options.rounds, options.dump);

  time_map["DeltaStepping"] =
      DeltaStepping_runner(G, options.src, options.rounds, options.ds_delta,
                           options.ds_num_buckets, options.dump);

  time_map["SSWidestPath"] = SSWidestPath_runner(
      G, options.src, options.rounds, options.wp_num_buckets,
      options.wp_no_blocked, options.wp_largemem, options.dump);

  for (const auto &[alg, time_per_iter] : time_map) {
    std::cout << "# # # " << alg << ", " << time_per_iter << "\n";
  }
}

} // namespace gbbs