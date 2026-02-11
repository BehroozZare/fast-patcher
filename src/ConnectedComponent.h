#pragma once
#include <vector>

// ---------------------------------------------------------------------------
// ConnectedComponent namespace: find connected components of a CSR graph.
//
// Graph convention:
//   Gp  -- row pointers   (size G_N + 1)
//   Gi  -- column indices  (size Gp[G_N])
//   G_N -- number of vertices
// ---------------------------------------------------------------------------
namespace ConnectedComponent {

// Compute connected components via BFS.
// Output: vertex_to_cc[v] = CC id for each vertex (sized to G_N).
// Returns the total number of connected components found.
int compute(
    const int* Gp, const int* Gi, int G_N,
    std::vector<int>& vertex_to_cc);

} // namespace ConnectedComponent
