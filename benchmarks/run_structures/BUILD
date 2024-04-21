licenses(["notice"])


cc_binary(
    name = "run_csr",
    srcs = ["run_csr.cc"],
    deps = ["//benchmarks:run_unweighted",
     ],
)

cc_binary(
    name = "run_gbbs",
    srcs = ["run_gbbs.cc"],
    deps = ["//benchmarks:run_unweighted",
    "//gbbs:bridge",
     ],
)

cc_binary(
    name = "run_csr_shuffled",
    srcs = ["run_csr.cc"],
    deps = ["//benchmarks:run_unweighted",
     ],
     defines = ["CSR_SHUFFLE"],
)

cc_binary(
    name = "run_std_set",
    srcs = ["run_std_set.cc"],
    deps = ["//benchmarks:run_unweighted",
     ],
)

cc_binary(
    name = "run_std_set_inplace",
    srcs = ["run_std_set.cc"],
    deps = ["//benchmarks:run_unweighted",
     ],
     defines = ["USE_INPLACE"],
)

cc_binary(
    name = "run_std_unordered_set",
    srcs = ["run_std_unordered_set.cc"],
    deps = ["//benchmarks:run_unweighted",
     ],
)

cc_binary(
    name = "run_std_unordered_set_inplace",
    srcs = ["run_std_unordered_set.cc"],
    deps = ["//benchmarks:run_unweighted",
     ],
     defines = ["USE_INPLACE"],
)

cc_binary(
    name = "run_absl_btree_set",
    srcs = ["run_absl_btree_set.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@absl//absl/container:btree",
     ],
)

cc_binary(
    name = "run_absl_btree_set_inplace",
    srcs = ["run_absl_btree_set.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@absl//absl/container:btree",
     ],
     defines = ["USE_INPLACE"],
)

cc_binary(
    name = "run_absl_flat_hash_set",
    srcs = ["run_absl_flat_hash_set.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@absl//absl/container:flat_hash_set",
     ],
)

cc_binary(
    name = "run_absl_flat_hash_set_inplace",
    srcs = ["run_absl_flat_hash_set.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@absl//absl/container:flat_hash_set",
     ],
     defines = ["USE_INPLACE"],
)

cc_binary(
    name = "run_vector_vector",
    srcs = ["run_vector_vector.cc"],
    deps = ["//benchmarks:run_unweighted",
    "//gbbs:bridge",
     ],
)

cc_binary(
    name = "run_sstgraph",
    srcs = ["run_sstgraph.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@SSTGraph//SSTGraph:SparseMatrix",
    "//gbbs:bridge",
     ],
)
cc_binary(
    name = "run_pcsr",
    srcs = ["run_pcsr.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@cpma//PMA:pcsr",
    "//gbbs:bridge",
     ],
)

cc_binary(
    name = "run_pcsr_orig",
    srcs = ["run_pcsr_orig.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@pcsr_orig//:pcsr_orig",
    "//gbbs:bridge",
     ],
)
cc_binary(
    name = "run_dhb",
    srcs = ["run_dhb.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@dhb//dhb:dynamic_hashed_blocks",
    "//gbbs:bridge",
    "@semisort//:semisort",
     ],
)

cc_binary(
    name = "run_terrace",
    srcs = ["run_terrace.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@Terrace//terrace:terrace_graph",
    "@ParallelTools//ParallelTools:parallel",
    "//gbbs:bridge",
     ],
)

cc_binary(
    name = "run_ppcsr",
    srcs = ["run_ppcsr.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@Terrace//terrace:PMA",
    "@ParallelTools//ParallelTools:parallel",
    "//gbbs:bridge",
    "@ParallelTools//ParallelTools:sort",
    "@ParallelTools//ParallelTools:reducer",
     ],
)

cc_binary(
    name = "run_vector_tinyset",
    srcs = ["run_vector_tinyset.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@SSTGraph//SSTGraph:TinySet",
    "//gbbs:bridge",
     ],
)

cc_binary(
    name = "run_vector_pma",
    srcs = ["run_vector_pma.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@cpma//PMA:cpma",
     ],
     defines = ["PMA_SETTINGS=pma_settings"],
)

cc_binary(
    name = "run_vector_spma",
    srcs = ["run_vector_pma.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@cpma//PMA:cpma",
     ],
     defines = ["PMA_SETTINGS=spmae_settings"],
)

cc_binary(
    name = "run_vector_cpma",
    srcs = ["run_vector_pma.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@cpma//PMA:cpma",
     ],
     defines = ["PMA_SETTINGS=cpma_settings"],
)

cc_binary(
    name = "run_vector_scpma",
    srcs = ["run_vector_pma.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@cpma//PMA:cpma",
     ],
     defines = ["PMA_SETTINGS=scpmae_settings"],
)

cc_binary(
    name = "run_single_pma",
    srcs = ["run_single_pma.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@cpma//PMA:cpma",
    "@reducer//:reducer",
    "//gbbs:bridge",
     ],
     defines = ["PMA_SETTINGS=pma_settings"],
)

cc_binary(
    name = "run_single_spma",
    srcs = ["run_single_pma.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@cpma//PMA:cpma",
    "@reducer//:reducer",
    "//gbbs:bridge",
     ],
     defines = ["PMA_SETTINGS=spmae_settings"],
)

cc_binary(
    name = "run_single_cpma",
    srcs = ["run_single_pma.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@cpma//PMA:cpma",
    "@reducer//:reducer",
    "//gbbs:bridge",
     ],
     defines = ["PMA_SETTINGS=cpma_settings"],
)

cc_binary(
    name = "run_single_scpma",
    srcs = ["run_single_pma.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@cpma//PMA:cpma",
    "@reducer//:reducer",
    "//gbbs:bridge",
     ],
     defines = ["PMA_SETTINGS=scpmae_settings"],
)

cc_binary(
    name = "run_vector_cpam",
    srcs = ["run_vector_cpam.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@CPAM//cpam",
    "//gbbs:bridge",
     ],
)

cc_binary(
    name = "run_vector_compressed_cpam",
    srcs = ["run_vector_cpam.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@CPAM//cpam",
    "//gbbs:bridge",
     ],
     defines = ["CPAM_COMPRESSED"],
)

cc_binary(
    name = "run_vector_compressed_cpam_inplace",
    srcs = ["run_vector_cpam.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@CPAM//cpam",
    "//gbbs:bridge",
     ],
     defines = ["CPAM_COMPRESSED", "USE_INPLACE"],
)

cc_binary(
    name = "run_vector_cpam_inplace",
    srcs = ["run_vector_cpam.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@CPAM//cpam",
    "//gbbs:bridge",
     ],
     defines = ["USE_INPLACE"],
)

cc_binary(
    name = "run_vector_aspen",
    srcs = ["run_vector_aspen.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@aspen//:aspen_lib",
    "//gbbs:macros",
     ],
)

cc_binary(
    name = "run_vector_compressed_aspen",
    srcs = ["run_vector_aspen.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@aspen//:aspen_lib",
    "//gbbs:macros",
     ],
    defines = ["ASPEN_USE_COMPRESSION"],
)

cc_binary(
    name = "run_vector_aspen_inplace",
    srcs = ["run_vector_aspen.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@aspen//:aspen_lib",
    "//gbbs:macros",
     ],
     defines = ["USE_INPLACE"],
)

cc_binary(
    name = "run_vector_compressed_aspen_inplace",
    srcs = ["run_vector_aspen.cc"],
    deps = ["//benchmarks:run_unweighted",
    "@aspen//:aspen_lib",
    "//gbbs:macros",
     ],
    defines = ["ASPEN_USE_COMPRESSION", "USE_INPLACE"],
)

package(
    default_visibility = ["//visibility:public"],
)
