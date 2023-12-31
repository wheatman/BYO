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

// Usage:
// numactl -i all ./BFS -src 10012 -s -m -rounds 3 twitter_SJ
// flags:
//   required:
//     -src: the source to compute the BFS from
//   optional:
//     -rounds : the number of times to run the algorithm
//     -c : indicate that the graph is compressed
//     -m : indicate that the graph should be mmap'd
//     -s : indicate that the graph is symmetric

#include "BFS.h"

namespace gbbs {

template <class Graph>
double BFS_runner(Graph& G, commandLine P) {
  uintE src = static_cast<uintE>(P.getOptionLongValue("-src", 0));
  std::cout << "### Application: BFS" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.N() << std::endl;
  std::cout << "### m: " << G.M() << std::endl;
  std::cout << "### Params: -src = " << src << std::endl;
  std::cout << "### ------------------------------------" << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  timer t;
  t.start();
  auto parents = BFS(G, src);
  double tt = t.stop();

  std::cout << "### Running Time: " << tt << std::endl;
#if 0
// useful for debugging
  std::vector<uintE> depths(G.N(), UINT_E_MAX);
  for (size_t j = 0; j < G.N(); j++) {
    uintE current_depth = 0;
    uintE current_parent = j;
    if (parents[j] == UINT_E_MAX) {
      continue;
    }
    while (current_parent != parents[current_parent]) {
      current_depth += 1;
      current_parent = parents[current_parent];
    }
    depths[j] = current_depth;
  }
  std::ofstream myfile;
  myfile.open("bfs.out");
  for (unsigned int i = 0; i < G.N(); i++) {
    myfile << depths[i] << std::endl;
  }
  myfile.close();
  #endif
  return tt;
}

}  // namespace gbbs

generate_main(gbbs::BFS_runner, false);
