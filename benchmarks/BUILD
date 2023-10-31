licenses(["notice"])

cc_library(
    name = "run_unweighted",
    hdrs = ["run_unweighted.h"],
    deps = ["//benchmarks/BFS/NonDeterministicBFS:BFS",
     "//benchmarks/ApproximateDensestSubgraph/ApproxPeelingBKV12:DensestSubgraph",
     "//benchmarks/ApproximateDensestSubgraph/GreedyCharikar:DensestSubgraph",
     "//benchmarks/CoSimRank:CoSimRank",
     "//benchmarks/Connectivity/BFSCC:Connectivity",
     "//benchmarks/Connectivity/LabelPropagation:Connectivity",
     "//benchmarks/Connectivity/SimpleUnionAsync:Connectivity",
     "//benchmarks/Connectivity/WorkEfficientSDB14:Connectivity",
     "//benchmarks/DegeneracyOrder/GoodrichPszona11:DegeneracyOrder",
     "//benchmarks/GraphColoring/Hasenplaugh14:GraphColoring",
     "//benchmarks/KCore/JulienneDBS17:KCore",
     "//benchmarks/LowDiameterDecomposition/MPX13:LowDiameterDecomposition",
     "//benchmarks/MaximalIndependentSet/RandomGreedy:MaximalIndependentSet",
     "//benchmarks/PageRank:PageRank",
     "//benchmarks/SSBetweenessCentrality/Brandes:SSBetweennessCentrality",
     "//benchmarks/Spanner/MPXV15:Spanner",
     "//benchmarks/Map/Base:Map",
     "//benchmarks/Map/RemoteRead:Map",
     ],
)


cc_binary(
    name = "run_unweighted_main",
    srcs = ["run_unweighted.cc"],
    deps = ["run_unweighted",
     ],
)
package(
    default_visibility = ["//visibility:public"],
)
