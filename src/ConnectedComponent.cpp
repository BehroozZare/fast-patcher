#include "ConnectedComponent.h"

#include <spdlog/spdlog.h>
// #include "BFS.h"  -- uncomment when implementing the TODO below

namespace ConnectedComponent {

// -----------------------------------------------------------------------
// Compute connected components via BFS.
//
// Walk every unvisited vertex, do a BFS from it, and label all reachable
// vertices with the same CC id.
// -----------------------------------------------------------------------
int compute(
    const int* Gp, const int* Gi, int G_N,
    std::vector<int>& vertex_to_cc)
{
    vertex_to_cc.clear();
    vertex_to_cc.resize(G_N, 0);
    spdlog::info("The number of connected component is 1, for G_N={}", vertex_to_cc.size());


    // TODO: implement BFS-based connected components using BFS::multi_source_bfs.
    //  - For each unvisited vertex v (vertex_to_cc[v] == -1):
    //      1. Set vertex_to_cc[v] = cc_id, prepare dist[v] = 0.
    //      2. Call BFS::multi_source_bfs(Gp, Gi, G_N, vertex_to_cc, dist)
    //         (single source -- only v has dist == 0).
    //         This will label all reachable vertices with cc_id.
    //      3. Increment cc_id.
    //  - Return the total number of CCs found.

    return 1;
}

} // namespace ConnectedComponent
