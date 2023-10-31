"""
This script runs all the benchmark binaries on specified input graphs.

Everything listed as a Bazel `cc_binary` under `benchmarks/` must be listed in
either `UNWEIGHTED_GRAPH_BENCHMARKS`, `WEIGHTED_GRAPH_BENCHMARKS`, or
`IGNORED_BINARIES` below. Forcing this explicit labeling stops contributors from
forgetting to update this file when they add a new benchmark.

The script only checks that the benchmarks run and exit without an error. It
does not check that the output of each benchmark is correct.

This script could be extended to further split benchmarks into ones that process
symmetric graphs versus asymmetric graphs, but currently the input graphs must
be symmetric so that all benchmarks can run on them.

This script should be invoked directly via Python >=3.7. Because this script
calls other Bazel commands, invoking it with `bazel run` won't work.
"""
from typing import List, Optional, Set, Tuple
import argparse
import fnmatch
import os
import sys
import subprocess

# The script will invoke these benchmark on an unweighted graph.
UNWEIGHTED_GRAPH_BENCHMARKS = [
    "//benchmarks/ApproximateDensestSubgraph/ApproxPeelingBKV12:DensestSubgraph_main",
    "//benchmarks/ApproximateDensestSubgraph/GreedyCharikar:DensestSubgraph_main",
    "//benchmarks/BFS/NonDeterministicBFS:BFS_main",
    "//benchmarks/CoSimRank:CoSimRank_main",
    "//benchmarks/Connectivity/BFSCC:Connectivity_main",
    "//benchmarks/Connectivity/LabelPropagation:Connectivity_main",
    "//benchmarks/Connectivity/SimpleUnionAsync:Connectivity_main",
    "//benchmarks/Connectivity/WorkEfficientSDB14:Connectivity_main",
    "//benchmarks/DegeneracyOrder/GoodrichPszona11:DegeneracyOrder_main",
    "//benchmarks/GraphColoring/Hasenplaugh14:GraphColoring_main",
    "//benchmarks/KCore/JulienneDBS17:KCore_main",
    "//benchmarks/LowDiameterDecomposition/MPX13:LowDiameterDecomposition_main",
    "//benchmarks/MaximalIndependentSet/RandomGreedy:MaximalIndependentSet_main",
    "//benchmarks/PageRank:PageRank_main",
    "//benchmarks/SSBetweenessCentrality/Brandes:SSBetweennessCentrality_main",
    "//benchmarks/Spanner/MPXV15:Spanner_main",
]

# The script will invoke these benchmarks on a weighted graph.
WEIGHTED_GRAPH_BENCHMARKS = [
    "//benchmarks/GeneralWeightSSSP/BellmanFord:BellmanFord_main",
    "//benchmarks/IntegralWeightSSSP/JulienneDBS17:wBFS_main",
    "//benchmarks/MinimumSpanningForest/Kruskal:MinimumSpanningForest_main",
    "//benchmarks/PositiveWeightSSSP/DeltaStepping:DeltaStepping_main",
    "//benchmarks/SSWidestPath/JulienneDBS17:SSWidestPath_main",
]

# The script will not invoke these binaries. Shell-style globbing is allowed in
# this list.
IGNORED_BINARIES = [
    # These incremental connectivity binaries take input in a format that's
    # different than that of other benchmarks.
    "//benchmarks/Connectivity/Incremental/mains:*_no_starting",
    "//benchmarks/SCAN/IndexBased/experiments:*",
    "//benchmarks/ApproximateDensestSubgraph/GreedyPlusPlus:DensestSubgraph_main",
    "//benchmarks/ApproximateSetCover/MANISBPT11:ApproximateSetCover_main",
    "//benchmarks/Biconnectivity/TarjanVishkin:Biconnectivity_main",
    "//benchmarks/CliqueCounting:Clique_main",
    "//benchmarks/Connectivity/ConnectIt/mains:unite_ldd",
    "//benchmarks/Connectivity/ConnectIt/mains:shiloach_vishkin",
    "//benchmarks/Connectivity/ConnectIt/mains:unite_rem_cas_bfs",
    "//benchmarks/Connectivity/ConnectIt/mains:unite_nd_nosample",
    "//benchmarks/Connectivity/ConnectIt/mains:unite_early_ldd",
    "//benchmarks/Connectivity/ConnectIt/mains:jayanti_nosample",
    "//benchmarks/Connectivity/ConnectIt/mains:unite_early_bfs",
    "//benchmarks/Connectivity/ConnectIt/mains:unite_rem_lock_nosample",
    "//benchmarks/Connectivity/ConnectIt/mains:liutarjan_nosample",
    "//benchmarks/Connectivity/ConnectIt/mains:unite_nd_bfs",
    "//benchmarks/Connectivity/ConnectIt/mains:unite_early_kout",
    "//benchmarks/Connectivity/ConnectIt/mains:unite_nosample",
    "//benchmarks/Connectivity/ConnectIt/mains:unite_nd_kout",
    "//benchmarks/Connectivity/ConnectIt/mains:unite_early_nosample",
    "//benchmarks/Connectivity/ConnectIt/mains:jayanti_kout",
    "//benchmarks/Connectivity/ConnectIt/mains:liutarjan_bfs",
    "//benchmarks/Connectivity/ConnectIt/mains:unite_rem_lock_kout",
    "//benchmarks/Connectivity/ConnectIt/mains:unite_bfs",
    "//benchmarks/Connectivity/ConnectIt/mains:liutarjan_kout",
    "//benchmarks/Connectivity/ConnectIt/mains:unite_rem_lock_ldd",
    "//benchmarks/Connectivity/ConnectIt/mains:unite_kout",
    "//benchmarks/DegeneracyOrder/BarenboimElkin08:DegeneracyOrder_main",
    "//benchmarks/Connectivity/Incremental/mains:jayanti_starting",
    "//benchmarks/Connectivity/Incremental/mains:liutarjan_starting",
    "//benchmarks/Connectivity/Incremental/mains:shiloachvishkin_starting",
    "//benchmarks/Connectivity/Incremental/mains:unite_early_starting",
    "//benchmarks/Connectivity/Incremental/mains:unite_nd_starting",
    "//benchmarks/Connectivity/Incremental/mains:unite_rem_cas_starting",
    "//benchmarks/Connectivity/Incremental/mains:unite_rem_lock_starting",
    "//benchmarks/Connectivity/Incremental/mains:unite_starting",
    "//benchmarks/Clustering/SeqHAC:HACSimilarity",
    "//benchmarks/Connectivity/ConnectIt/mains:jayanti_ldd",
    "//benchmarks/Connectivity/ConnectIt/mains:bfscc",
    "//benchmarks/Clustering/SeqHAC:HACDissimilarity",
    "//benchmarks/Connectivity/ConnectIt/mains:unite_rem_cas_nosample",
    "//benchmarks/Connectivity/ConnectIt/mains:jayanti_bfs",
    "//benchmarks/Connectivity/ConnectIt/mains:unite_rem_cas_kout",
    "//benchmarks/Connectivity/ConnectIt/mains:unite_rem_cas_ldd",
    "//benchmarks/Connectivity/ConnectIt/mains:gbbscc",
    "//benchmarks/Connectivity/ConnectIt/mains:liutarjan_ldd",
    "//benchmarks/Connectivity/ConnectIt/mains:unite_nd_ldd",
    "//benchmarks/Connectivity/ConnectIt/mains:unite_rem_lock_bfs",
    "//benchmarks/Connectivity/ConnectIt/mains:label_propagation",
    "//benchmarks/CycleCounting/Parallel5Cycle:FiveCycle_main",
    "//benchmarks/KCore/ApproximateKCore:KCore_main",
    "//benchmarks/KTruss:KTruss_main",
    "//benchmarks/MaximalIndependentSet/Yoshida:MaximalIndependentSet_main",
    "//benchmarks/MaximalMatching/RandomGreedy:MaximalMatching_main",
    "//benchmarks/MaximalMatching/Yoshida:MaximalMatching_main",
    "//benchmarks/SCAN/IndexBased:SCAN_main",
    "//benchmarks/TriangleCounting/ShunTangwongsan15:Triangle_main",
    "//benchmarks/SpanningForest/LabelPropagation:SpanningForest_main",
    "//benchmarks/SpanningForest/BFSSF:SpanningForest_main",
    "//benchmarks/Connectivity/ConnectIt/mains:label_propagation",
    "//benchmarks/SpanningForest/SDB14:SpanningForest_main",
    "//benchmarks/MinimumSpanningForest/Boruvka:MinimumSpanningForest_main",
    "//benchmarks/StronglyConnectedComponents/RandomGreedyBGSS16:StronglyConnectedComponents_main",
]


def get_all_benchmark_binaries() -> List[str]:
    """Returns a list of all binaries under `benchmarks/`."""
    return subprocess.run(
        ["bazel", "query", "kind(cc_binary, //benchmarks/...)"],
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    ).stdout.splitlines()


def check_listed_binaries(
    valid_binaries: List[str], ignored_binaries: List[str]
) -> None:
    """Checks listed binaries for consistency.

    Args:
        valid_binaries: Names of binaries that are considered valid.
        ignored_binaries: Names of binaries that are considered invalid
            and should not be run.

    Raises:
        ValueError: `valid_binaries` and `ignored_binaries` overlap.
        ValueError: `valid_binaries` and `ignored_binaries` combined don't equal
            the set of all benchmark-related binaries as determined by
            `get_all_benchmark_binaries()`.
    """
    remaining_binaries = set(get_all_benchmark_binaries())
    for ignored_binaries_pattern in ignored_binaries:
        conflicting_binaries = fnmatch.filter(
            valid_binaries, ignored_binaries_pattern)
        if conflicting_binaries:
            raise ValueError(
                "Benchmarks listed as both valid and ignored: {}".format(
                    conflicting_binaries
                )
            )
        binaries_to_ignore = fnmatch.filter(
            remaining_binaries, ignored_binaries_pattern
        )
        if not binaries_to_ignore:
            print(
                "Warning: ignore rule {} has no effect".format(
                    ignored_binaries_pattern)
            )
        remaining_binaries -= set(binaries_to_ignore)

    valid_binaries = set(valid_binaries)
    if remaining_binaries != valid_binaries:
        extra_listed_binaries = valid_binaries - remaining_binaries
        if extra_listed_binaries:
            raise ValueError(
                "Listed benchmarks do not exist: {}".format(
                    extra_listed_binaries)
            )
        missing_listed_binaries = remaining_binaries - valid_binaries
        if missing_listed_binaries:
            raise ValueError(
                "Please update {} to include binaries {}".format(
                    __file__, missing_listed_binaries
                )
            )


def run_all_benchmarks(
    graph_benchmarks: List[str],
    build_options: str,
    graph_file: str,
    are_graphs_compressed: bool,
    timeout: int,
    weighted_graph: bool,
    src_vertex: Optional[str]
) -> List[Tuple[str, str]]:
    """Runs all benchmarks, returning a list of failing benchmarks.

    Args:
        graph_benchmarks: List of all benchmarks to run on the
            graph.
        build_options: list of options to bass through the the buod with --copt=build_options
        graph_file: File path to the  graph.
        are_compressed_compressed: Whether the graph files hold compressed
            graphs.
        timeout: Benchmarks that run longer than this timeout period in seconds
            are considered to have failed. If this is `None` then the benchmarks
            have no time limit.
        weighted_graph: if the graph is weighted
        src_vertex what source vrtex to use for benchmarks that require a source

    Returns:
        A list of names of benchmarks that fail along with a failure reason.
    """

    BAZEL_FLAGS = ["--compilation_mode", "opt"]
    gbbs_flags = ["-s", "-rounds", "10", "-src", src_vertex]
    if are_graphs_compressed:
        gbbs_flags += ["-c"]

    benchmarks = graph_benchmarks
    # Compile all the benchmarks up front --- it's faster than compiling
    # them individually since Bazel can compile several files in parallel.
    command = " ".join(["bazel", "build"] + BAZEL_FLAGS +
                   ["--keep_going", "--output_filter=DONT_MATCH_ANYTHING"] +  [build_options] + benchmarks)
    print(command)

    subprocess.run(command, shell=True)

    failed_benchmarks = []

    def test_benchmark(
        benchmark: str, graph_file: str, additional_gbbs_flags: List[str]
    ) -> None:
        try:
            command = " ".join(["bazel", "run"]
                + BAZEL_FLAGS
                + [benchmark, "--"]
                + gbbs_flags
                + additional_gbbs_flags
                + [graph_file])
            benchmark_run = subprocess.run(command, shell=True,
                timeout=timeout, capture_output=True,
            )
            if benchmark_run.returncode:
                print(benchmark_run)
                failed_benchmarks.append(
                    (
                        benchmark,
                        "Exited with error code {}".format(
                            benchmark_run.returncode),
                    )
                )
            else:
                print(benchmark, benchmark_run.stdout.decode("utf-8").split("\n")[-2].split(":")[-1])
        except subprocess.TimeoutExpired:
            failed_benchmarks.append((benchmark, "Timeout"))

    if weighted_graph:
        extra_gbbs_flags = ["-w", build_options]
    else:
        extra_gbbs_flags = [build_options]
    
    for benchmark in graph_benchmarks:
        test_benchmark(
            benchmark=benchmark,
            graph_file=graph_file,
            additional_gbbs_flags = extra_gbbs_flags
        )

    return failed_benchmarks

def get_build_options():
    build_options = []
    build_options.append("")

    # # vary implementations

    build_options.append("--copt='-DUNWEIGHTED_SYM_GRAPH_IMPL=gbbs::graph_implementations::symmetric_graph<gbbs::symmetric_vertex,gbbs::empty,true>'")
    build_options.append("--copt='-DUNWEIGHTED_SYM_GRAPH_IMPL=gbbs::graph_implementations::symmetric_set_graph<gbbs::symmetric_vertex,gbbs::empty,gbbs::graph_implementations::vector_set>'")

    build_options.append("--copt='-DUNWEIGHTED_SYM_GRAPH_IMPL=gbbs::graph_implementations::symmetric_set_graph<gbbs::symmetric_vertex,gbbs::empty,std::set<gbbs::uintE>>'")
    build_options.append("--copt='-DUNWEIGHTED_SYM_GRAPH_IMPL=gbbs::graph_implementations::symmetric_set_graph<gbbs::symmetric_vertex,gbbs::empty,std::unordered_set<gbbs::uintE>>'")
    build_options.append("--copt='-DUNWEIGHTED_SYM_GRAPH_IMPL=gbbs::graph_implementations::symmetric_set_graph<gbbs::symmetric_vertex,gbbs::empty,absl::btree_set<gbbs::uintE>>'")
    build_options.append("--copt='-DUNWEIGHTED_SYM_GRAPH_IMPL=gbbs::graph_implementations::symmetric_set_graph<gbbs::symmetric_vertex,gbbs::empty,absl::flat_hash_set<gbbs::uintE>>'")

    build_options.append("--copt='-DUNWEIGHTED_SYM_GRAPH_IMPL=gbbs::graph_implementations::symmetric_SSTGraph_graph<gbbs::symmetric_vertex,gbbs::empty>'")

    build_options.append("--copt='-DUNWEIGHTED_SYM_GRAPH_IMPL=gbbs::graph_implementations::symmetric_dhb_graph<gbbs::symmetric_vertex,gbbs::empty>'")

    # vary API

    build_options.append("--copt='-DGRAPH_API_DEGREE=false' --copt='-DGRAPH_API_M=false' --copt='-DGRAPH_API_PARALLEL_MAP=false' --copt='-DGRAPH_API_MAP_EARLY_EXIT=false' --copt='-DGRAPH_API_PARALLEL_MAP_EARLY_EXIT=false'")

    build_options.append("--copt='-DGRAPH_API_DEGREE=true' --copt='-DGRAPH_API_M=true' --copt='-DGRAPH_API_PARALLEL_MAP=true' --copt='-DGRAPH_API_MAP_EARLY_EXIT=false' --copt='-DGRAPH_API_PARALLEL_MAP_EARLY_EXIT=false'")

    build_options.append("--copt='-DGRAPH_API_DEGREE=true' --copt='-DGRAPH_API_M=true' --copt='-DGRAPH_API_PARALLEL_MAP=true' --copt='-DGRAPH_API_MAP_EARLY_EXIT=true' --copt='-DGRAPH_API_PARALLEL_MAP_EARLY_EXIT=false'")

    build_options.append("--copt='-DGRAPH_API_DEGREE=true' --copt='-DGRAPH_API_M=true' --copt='-DGRAPH_API_PARALLEL_MAP=false' --copt='-DGRAPH_API_MAP_EARLY_EXIT=true' --copt='-DGRAPH_API_PARALLEL_MAP_EARLY_EXIT=false'")

    build_options.append("--copt='-DGRAPH_API_DEGREE=false' --copt='-DGRAPH_API_M=true' --copt='-DGRAPH_API_PARALLEL_MAP=true' --copt='-DGRAPH_API_MAP_EARLY_EXIT=true' --copt='-DGRAPH_API_PARALLEL_MAP_EARLY_EXIT=true'")

    build_options.append("--copt='-DGRAPH_API_DEGREE=true' --copt='-DGRAPH_API_M=false' --copt='-DGRAPH_API_PARALLEL_MAP=true' --copt='-DGRAPH_API_MAP_EARLY_EXIT=true' --copt='-DGRAPH_API_PARALLEL_MAP_EARLY_EXIT=true'")

    build_options.append("--copt='-DGRAPH_API_DEGREE=true' --copt='-DGRAPH_API_M=false' --copt='-DGRAPH_API_PARALLEL_MAP=false' --copt='-DGRAPH_API_MAP_EARLY_EXIT=false' --copt='-DGRAPH_API_PARALLEL_MAP_EARLY_EXIT=false'")

    return build_options


if __name__ == "__main__":
    check_listed_binaries(
        valid_binaries=UNWEIGHTED_GRAPH_BENCHMARKS + WEIGHTED_GRAPH_BENCHMARKS,
        ignored_binaries=IGNORED_BINARIES,
    )

    parser = argparse.ArgumentParser(
        description=(
            "Runs all benchmarks on the specified input graphs to check "
            "whether the benchmarks run and exit without errors."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--graph",
        "-g",
        type=str,
        help=(
            "Absolute path to an graph on which to run all "
            "graph benchmarks. If not provided, those benchmarks "
            "will not be run."
        ),
    )
    parser.add_argument(
        "--weighted",
        "-w",
        action="store_true",
        help="Add this flag if input graphs are weighted.",
    )
    parser.add_argument(
        "--compressed",
        "-c",
        action="store_true",
        help="Add this flag if input graphs are compressed graphs.",
    )
    parser.add_argument(
        "--timeout",
        "-t",
        type=float,
        default=600,
        help="(seconds) - Halt benchmarks that run longer than this time.",
    )
    parser.add_argument(
        "--src_vertex",
        type=str,
        default="0",
        help="which source vertex to use for algorithms that use a source",
    )
    parsed_args = parser.parse_args()


    graph_file = (
        os.path.abspath(parsed_args.graph)
    )

    benchmarks = UNWEIGHTED_GRAPH_BENCHMARKS if not parsed_args.weighted else WEIGHTED_GRAPH_BENCHMARKS

    for build_option in get_build_options():
        print(build_option)
        failed_benchmarks = run_all_benchmarks(
            graph_benchmarks=benchmarks,
            build_options = build_option,
            graph_file=graph_file,
            are_graphs_compressed=parsed_args.compressed,
            timeout=parsed_args.timeout,
            weighted_graph = parsed_args.weighted,
            src_vertex=parsed_args.src_vertex
        )
        if failed_benchmarks:
            print("Benchmarks failed: {}".format(failed_benchmarks))
        else:
            print("Success! All benchmarks completed without an error.")
