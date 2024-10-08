# BYO: A Unified Framework for Benchmarking Large-Scale Graph Containers

Organization
--------

This repository contains code for our paper "[BYO: A Unified Framework for Benchmarking Large-Scale Graph Containers](https://arxiv.org/abs/2405.11671)" (VLDB'24).

It is designed to make it as easy as possible to implement and benchmark new graph container data structures.

The containers already implemented can be found in benchmarks/run_structures and include:

The following set based structures:
- vector of absl::btree [code](https://github.com/abseil/abseil-cpp/tree/master)
  - run_absl_btree_set
  - run_absl_btree_set_inplace
- vector of absl::flat_hash_map [code](https://github.com/abseil/abseil-cpp/tree/master)
  - run_absl_flat_hash_set
  - run_absl_flat_hash_set_inplace
- vector of std::set
  - run_std_set
  - run_std_set_inplace
- vector of std::unordered_set
  - run_std_unordered_set
  - run_std_unordered_set
- vector of aspen trees [code](https://github.com/DapperX/aspen/) [paper](https://dl.acm.org/doi/10.1145/3314221.3314598)
  - run_vector_aspen
  - run_vector_compressed_aspen (as described in the paper)
  - run_vector_aspen_inplace
  - run_vector_compressed_aspen_inplace
- vector of cpam trees [code](https://github.com/ParAlg/CPAM) [paper](https://dl.acm.org/doi/abs/10.1145/3519939.3523733)
  - run_vector_cpam
  - run_vector_compressed_cpam (as described in the paper)
  - run_vector_cpam_inplace
  - run_vector_compressed_cpam_inplace
- vector of PMAs [code](https://github.com/wheatman/Packed-Memory-Array/)
  - run_vector_pma
  - run_vector_spma
  - run_vector_cpma
  - run_vector_scpma
- vector of tinysets [code](https://github.com/wheatman/SSTGraph/tree/main) [paper](https://ieeexplore.ieee.org/abstract/document/9671836)
  - run_vector_tinyset
- vector of vectors
  - run_vector_vector


The following full graph based structures:
- csr
  - run_csr
  - run_csr_shuffled
- gbbs [paper](https://dl.acm.org/doi/10.1145/3434393)
  - run_gbbs (as described in the paper)
- dhb [code](https://github.com/wheatman/dhb/tree/main) [paper](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.SEA.2022.11)
  - run_dhb (as described in the paper)
- pcsr original [code](https://github.com/wheatman/Packed-Compressed-Sparse-Row) [paper](https://ieeexplore.ieee.org/abstract/document/8547566)
  - run_pcsr_orig (can only run small graphs)
- pcsr [code](https://github.com/wheatman/Packed-Memory-Array/)
  - run_pcsr (as described in the paper)
- ppcsr [code](https://github.com/PASSIONLab/terrace) [paper](https://epubs.siam.org/doi/abs/10.1137/1.9781611976472.3)
  - run_ppcsr (as described in the paper)
- F-Graph [code](https://github.com/wheatman/Packed-Memory-Array/) [paper](https://dl.acm.org/doi/abs/10.1145/3627535.3638492)
  - run_single_pma
  - run_single_spma
  - run_single_cpma 
  - run_single_scpma (as described in the paper)
- SSTGraph [code](https://github.com/wheatman/SSTGraph) [paper](https://ieeexplore.ieee.org/abstract/document/9671836)
  - run_sstgraph (as described in the paper)
- Terrace [code](https://github.com/PASSIONLab/terrace) [paper](https://dl.acm.org/doi/abs/10.1145/3448016.3457313)
  - run_terrace (as described in the paper)

For many of the data structures we implement multiple different versions, we note above which version was described and benchmarked in the papers.

The raw data for all systems we tested can be found at [here](https://docs.google.com/spreadsheets/d/1Vi3bbCeWBCgl-Me15aAYCf0wPGSVTteuVJFeKkLjcW8/edit?usp=sharing)


Compilation
--------

Compiler:
* g++ &gt;= 11 with pthread support (Homemade Scheduler)
* clang++ &gt;= 14 with support for OpenCilk

Build system:
* [Bazel](https://docs.bazel.build/versions/master/install.html) 2.1.0


The default compilation uses a lightweight scheduler developed at CMU (Homemade)
for parallelism.


To compile codes for graphs with more than 2^32 edges, the `LONG` command-line
parameter should be set. If the graph has more than 2^32 vertices, the
`EDGELONG` command-line parameter should be set. Note that the codes have not
been tested with more than 2^32 vertices, so if any issues arise please contact
make an issue on github.

To compile with the OpenCilk scheduler instead of the Homegrown scheduler, use
the Bazel configuration `--config=cilk`. To compile using OpenMP instead, use
the Bazel configuration `--config=openmp`. To compile serially instead, use the
Bazel configuration `--config=serial`. 

To build:
```sh
# Load external libraries as submodules. (This only needs to be run once.)
git submodule update --init

# Note that the default compilation mode in bazel is to build optimized binaries
# (stripped of debug symbols). You can compile debug binaries by supplying `-c
# dbg` to the bazel build command.

#To build all of the different structures use the command 
`bazel build benchmarks/run_structures:all`

#similarly individual systems can be built with 

bazel build benchmarks/run_structures:run_csr

# The following commands cleans the directory:

$ bazel clean  # removes all executables

```


Adding a new graph Container
-------
To add a new set container we recommend following the example given in run_vector.

Similarly to add a new Graph container we recommend following the example in run_single_pma.

You will also need to add the new container to the build file in `benchmarks/run_structures/BUILD` so that it can be built.



Running code
-------
The applications take the input graph as input as well as an optional
flag "-s" to indicate a symmetric graph.  Symmetric graphs should be
called with the "-s" flag for better performance. For example:

```sh
# For Bazel:
$ ./bazel-bin/benchmarks/run_structures/run_csr -s -src 10 ~/gbbs/inputs/rMatGraph_J_5_100
```

Note that the codes that compute single-source shortest paths (or centrality)
take an extra `-src` flag. The benchmark is run four times by default, and can
be changed by passing the `-rounds` flag followed by an integer indicating the
number of runs.

On NUMA machines, adding the command "numactl -i all " when running
the program may improve performance for large graphs. For example:

```sh
$ numactl -i all ./bazel-bin/benchmarks/run_structures/run_csr [...]
```

Running code on compressed graphs
-----------
The CSR graphs can also take in a compressed graph.  The dynamic systems take in an uncompressed graph and compress it on the fly

We make use of the bytePDA format in our benchmark, which is similar to the
parallelByte format of Ligra+, extended with additional functionality. We have
provided a converter utility which takes as input an uncompressed graph and
outputs a bytePDA graph. The converter can be used as follows:

```sh
# For Bazel:
bazel run //utils:compressor -- -s -o ~/gbbs/inputs/rMatGraph_J_5_100.bytepda ~/gbbs/inputs/rMatGraph_J_5_100
bazel run //utils:compressor -- -s -w -o ~/gbbs/inputs/rMatGraph_WJ_5_100.bytepda ~/gbbs/inputs/rMatGraph_WJ_5_100

```

After an uncompressed graph has been converted to the bytepda format,
applications can be run on it by passing in the usual command-line flags, with
an additional `-c` flag.

```sh
# For Bazel:
$ bazel run //benchmarks/BFS/NonDeterministicBFS:BFS_main -- -s -c -src 10 ~/gbbs/inputs/rMatGraph_J_5_100.bytepda

```

When processing large compressed graphs, using the `-m` command-line flag can
help if the file is already in the page cache, since the compressed graph data
can be mmap'd. Application performance will be affected if the file is not
already in the page-cache. We have found that using `-m` when the compressed
graph is backed by SSD results in a slow first-run, followed by fast subsequent
runs.

Running code on binary-encoded graphs
-----------
We make use of a binary-graph format in our benchmark. The binary representation
stores the representation we use for in-memory processing (compressed sparse row)
directly on disk, which enables applications to avoid string-conversion overheads
associated with the adjacency graph format described below. We have provided a
converter utility which takes as input an uncompressed graph (e.g., in adjacency
graph format) and outputs this graph in the binary format. The converter can be
used as follows:

```sh
# For Bazel:
bazel run //utils:compressor -- -s -o ~/gbbs/inputs/rMatGraph_J_5_100.binary ~/gbbs/inputs/rMatGraph_J_5_100

```

After an uncompressed graph has been converted to the binary format,
applications can be run on it by passing in the usual command-line flags, with
an additional `-b` flag. Note that the application will always load the binary
file using mmap.

```sh
# For Bazel:
$ bazel run //benchmarks/BFS/NonDeterministicBFS:BFS_main -- -s -b -src 10 ~/gbbs/inputs/rMatGraph_J_5_100.binary

```

Note that application performance will be affected if the file is not already
in the page-cache. We have found that using `-m` when the binary graph is backed
by SSD or disk results in a slow first-run, followed by fast subsequent runs.


Input Formats
-----------
We support the adjacency graph format used by the [Problem Based Benchmark
suite](http://www.cs.cmu.edu/~pbbs/benchmarks/graphIO.html)
and [Ligra](https://github.com/jshun/ligra).

The adjacency graph format starts with a sequence of offsets one for each
vertex, followed by a sequence of directed edges ordered by their source vertex.
The offset for a vertex i refers to the location of the start of a contiguous
block of out edges for vertex i in the sequence of edges. The block continues
until the offset of the next vertex, or the end if i is the last vertex. All
vertices and offsets are 0 based and represented in decimal. The specific format
is as follows:

```
AdjacencyGraph
<n>
<m>
<o0>
<o1>
...
<o(n-1)>
<e0>
<e1>
...
<e(m-1)>
```

This file is represented as plain text.

Weighted graphs are represented in the weighted adjacency graph format. The file
should start with the string "WeightedAdjacencyGraph". The m edge weights
should be stored after all of the edge targets in the .adj file.

**Using SNAP graphs**

Graphs from the [SNAP dataset
collection](https://snap.stanford.edu/data/index.html) are commonly used for
graph algorithm benchmarks. We provide a tool that converts the most common SNAP
graph format to the adjacency graph format that GBBS accepts. Usage example:
```sh
# Download a graph from the SNAP collection.
wget https://snap.stanford.edu/data/wiki-Vote.txt.gz
gzip --decompress ${PWD}/wiki-Vote.txt.gz
# Run the SNAP-to-adjacency-graph converter.
# Run with Bazel:
bazel run //utils:snap_converter -- -s -i ${PWD}/wiki-Vote.txt -o <output file>
# Or run with Make:
#   cd utils
#   make snap_converter
#   ./snap_converter -s -i <input file> -o <output file>
```

## Citation
Please cite as 

Brian Wheatman, Xiaojun Dong, Zheqi Shen, Laxman Dhulipala, Jakub Łącki, Prashant Pandey, and Helen Xu. 2024. BYO: A Unified Framework for Benchmarking Large-Scale Graph Containers. Proc. VLDB Endow. 17, 9 (May 2024), 2307–2320. https://doi.org/10.14778/3665844.3665859
