#include "graph_io.h"

#include <sys/mman.h>

namespace gbbs {
namespace gbbs_io {

typedef std::pair<uintE, uintE> intPair;

namespace internal {

// Starting from the current position, skips all consecutive lines in the stream
// that start with '#' or are empty.
//
// The intent here is that lines starting with '#' are interpreted to be
// comments that should be ignored.
void skip_ifstream_comments(std::ifstream* stream) {
  std::string line;
  while (*stream) {
    std::streampos current_position = stream->tellg();
    std::getline(*stream, line);
    if (!(line.empty() || line[0] == '#')) {
      stream->seekg(current_position);
      return;
    }
  }
}

}  // namespace internal

template <>
Edge<gbbs::empty>::Edge(const uintE _from, const uintE _to)
    : from(_from), to(_to) {}

std::tuple<size_t, size_t, uintT*, uintE*> parse_unweighted_graph(
    const char* fname, bool mmap, bool binary, char* bytes, size_t bytes_size) {
  uintT* offsets;
  uintE* edges;
  uint64_t n, m;

  if (!binary) {
    sequence<char> S;

    if (bytes == nullptr) {
      if (mmap) {
        std::pair<char*, size_t> MM = mmapStringFromFile(fname);
        S = sequence<char>(MM.second);
        // Cannot mutate the graph unless we copy.
        parallel_for(0, S.size(), [&](size_t i) { S[i] = MM.first[i]; });
        if (munmap(MM.first, MM.second) == -1) {
          perror("munmap");
          exit(-1);
        }
      } else {
        S = readStringFromFile(fname);
      }
    }
    sequence<slice<char>> tokens = parlay::map_tokens(
        parlay::make_slice(S), [](auto x) { return parlay::make_slice(x); });

    debug(std::string header = std::string(tokens[0].begin(), tokens[0].size());
          assert(header == internal::kUnweightedAdjGraphHeader););

    n = parlay::internal::chars_to_int_t<unsigned long>(tokens[1]);
    m = parlay::internal::chars_to_int_t<unsigned long>(tokens[2]);

    debug(std::cout << "# n = " << n << " m = " << m
                    << " len = " << (tokens.size() - 1) << "\n";
          uint64_t len = tokens.size() - 1; assert(len == n + m + 2););

    offsets = gbbs::new_array_no_init<uintT>(n + 1);
    edges = gbbs::new_array_no_init<uintE>(m);

    parallel_for(0, n, [&](size_t i) {
      offsets[i] =
          parlay::internal::chars_to_int_t<unsigned long>(tokens[i + 3]);
    });
    offsets[n] = m; /* make sure to set the last offset */
    parallel_for(0, m, [&](size_t i) {
      edges[i] =
          parlay::internal::chars_to_int_t<unsigned long>(tokens[i + n + 3]);
    });

    S.clear();
    tokens.clear();
  } else {
    std::pair<char*, size_t> MM = mmapStringFromFile(fname);

    auto mmap_file = MM.first;

    long* sizes = (long*)mmap_file;
    n = sizes[0], m = sizes[1];

    offsets = (uintT*)(mmap_file + 3 * sizeof(long));
    uint64_t skip = 3 * sizeof(long) + (n + 1) * sizeof(intT);
    edges = (uintE*)(mmap_file + skip);
  }

  return std::make_tuple(n, m, offsets, edges);
}


std::tuple<char*, size_t> parse_compressed_graph(const char* fname, bool mmap) {
  char* bytes;
  size_t bytes_size;

  if (mmap) {
    std::tie(bytes, bytes_size) = mmapStringFromFile(fname);
  } else {
    std::tie(bytes, bytes_size) = read_o_direct(fname);
  }
  return std::make_tuple(bytes, bytes_size);
}

std::vector<Edge<gbbs::empty>> read_unweighted_edge_list(const char* filename) {
  std::ifstream file{filename};
  if (!file.is_open()) {
    std::cout << "ERROR: Unable to open file: " << filename << '\n';
    std::terminate();
  }
  internal::skip_ifstream_comments(&file);

  std::vector<Edge<gbbs::empty>> edge_list;
  uintE from;
  uintE to;
  while (file >> from >> to) {
    edge_list.emplace_back(from, to);
  }
  return edge_list;
}

}  // namespace gbbs_io
}  // namespace gbbs
