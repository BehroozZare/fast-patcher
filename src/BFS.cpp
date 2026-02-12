#include "BFS.h"
// When implementing the TODOs below you will likely need:
//   #include <queue>

namespace BFS {

// -----------------------------------------------------------------------
// Multi-source BFS on the full graph.
//
// Precondition: caller has set label[s] and dist[s] = 0 for each source s.
//               All other entries should be -1.
// -----------------------------------------------------------------------
void multi_source_bfs(
    const int* Gp, const int* Gi, int G_N,
    const std::vector<int>& source_vertices,
    std::vector<int>& label,
    std::vector<int>& dist)
{
    // TODO: implement multi-source BFS.
    //  - Enqueue every vertex v where dist[v] == 0.
    //  - BFS level by level.  For each dequeued vertex v, iterate over
    //    its neighbours u (Gi[Gp[v]] .. Gi[Gp[v+1]-1]):
    //      * If label[u] == -1, set label[u] = label[v],
    //        dist[u] = dist[v] + 1, and enqueue u.
}

// -----------------------------------------------------------------------
// Restricted multi-source BFS.
//
// Same as multi_source_bfs but only expands into vertices where
// restrict_mask[v] == cluster_id.  Vertices outside the cluster are
// never enqueued.
// -----------------------------------------------------------------------
void restricted_bfs(
    const int* Gp, const int* Gi, int G_N,
    const std::vector<int>& restrict_mask,
    int cluster_id,
    std::vector<int>& label,
    std::vector<int>& dist)
{
    // TODO: implement restricted multi-source BFS.
    //  - Enqueue every vertex v where dist[v] == 0
    //    (and restrict_mask[v] == cluster_id).
    //  - BFS level by level.  For each dequeued vertex v, iterate over
    //    its neighbours u:
    //      * Skip u if restrict_mask[u] != cluster_id.
    //      * If label[u] == -1, set label[u] = label[v],
    //        dist[u] = dist[v] + 1, and enqueue u.
}

} // namespace BFS
