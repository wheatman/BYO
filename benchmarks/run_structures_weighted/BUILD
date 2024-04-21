licenses(["notice"])


cc_binary(
    name = "run_csr",
    srcs = ["run_csr.cc"],
    deps = ["//benchmarks:run_weighted",
     ],
)

cc_binary(
    name = "run_std_map",
    srcs = ["run_simple_map.cc"],
    deps = ["//benchmarks:run_weighted",
     ],
     defines = ["STD_MAP"],
)

cc_binary(
    name = "run_std_unordered_map",
    srcs = ["run_simple_map.cc"],
    deps = ["//benchmarks:run_weighted",
     ],
     defines = ["STD_UNORDERED_MAP"],
)

cc_binary(
    name = "run_absl_btree_map",
    srcs = ["run_simple_map.cc"],
    deps = ["//benchmarks:run_weighted",
            "@absl//absl/container:btree",
     ],
     defines = ["ABSL_BTREE"],
)

cc_binary(
    name = "run_absl_flat_map",
    srcs = ["run_simple_map.cc"],
    deps = ["//benchmarks:run_weighted",
            "@absl//absl/container:flat_hash_map",
     ],
     defines = ["ABSL_FLAT"],
)

package(
    default_visibility = ["//visibility:public"],
)
