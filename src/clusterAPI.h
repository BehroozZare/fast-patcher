#pragma once
#include <vector>
#include "Lloyd.h"

namespace Patcher {

// Main entry point: partition graph into clusters of size <= patch_size.
//
// Graph is given as raw CSR arrays:
//   Gp  -- row pointers   (size n + 1)
//   Gi  -- column indices  (size Gp[n])
//   n   -- number of vertices
//
// opt can be nullptr (defaults will be used).
//
// Output: vertex_to_cluster[v] = cluster id for every vertex v.
void create_clusters(
    int n,
    const int* Gp,
    const int* Gi,
    int patch_size,
    const LloydOptions* opt,
    std::vector<int>& vertex_to_cluster);

} // namespace Patcher
