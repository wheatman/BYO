cc_library(
  name = "KCore",
  hdrs = ["KCore.h"],
  deps = [
  "//gbbs:gbbs",
  "//gbbs:julienne",
  ]
)

cc_binary(
  name = "KCore_main",
  srcs = ["KCore.cc"],
  deps = [
  ":KCore",
  "//benchmarks/KCore/JulienneDBS17:KCore",
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
