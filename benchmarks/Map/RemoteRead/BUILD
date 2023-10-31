licenses(["notice"])

cc_library(
    name = "Map",
    hdrs = ["Map.h"],
    deps = ["//gbbs",
            "@reducer//:reducer"
            ],
)

cc_binary(
    name = "Map_main",
    srcs = ["Map.cc"],
    deps = [":Map"],
)

package(
    default_visibility = ["//visibility:public"],
)
