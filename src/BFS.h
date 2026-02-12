#pragma once
#include <vector>

// ---------------------------------------------------------------------------
// BFS namespace: reusable BFS primitives on CSR graphs.
//
// Graph convention:
//   Gp  -- row pointers   (size G_N + 1)
//   Gi  -- column indices  (size Gp[G_N])
//   G_N -- number of vertices
//
// Both functions expect the caller to set up label[] and dist[] before
// calling.  Specifically:
//   - For each source vertex s: label[s] = <source id>, dist[s] = 0.
//   - For all other vertices:   label[v] = -1,           dist[v] = -1.
// ---------------------------------------------------------------------------
namespace BFS {

// Multi-source BFS on the full graph.
// Expands from every vertex where dist[v] == 0, level by level.
// After return:
//   label[v] = label of the source that first reached v.
//   dist[v]  = BFS distance from v to that source.
void multi_source_bfs(
    const int* Gp, const int* Gi, int G_N,
    const std::vector<int>& source_vertices,
    std::vector<int>& label,
    std::vector<int>& dist);

// Restricted multi-source BFS.
// Same as above, but only expands into vertices where
// restrict_mask[v] == cluster_id.  All other vertices are skipped.
// Useful for intra-cluster BFS (e.g., from boundary inward).
void restricted_bfs(
    const int* Gp, const int* Gi, int G_N,
    const std::vector<int>& restrict_mask,
    int cluster_id,
    std::vector<int>& label,
    std::vector<int>& dist);

} // namespace BFS
