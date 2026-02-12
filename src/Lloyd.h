#pragma once
#include <vector>

// ---------------------------------------------------------------------------
// Configuration for the Lloyd clustering.
// Matches the "opt" struct described in ALGORITHM_STEPS.md.
// ---------------------------------------------------------------------------
struct LloydOptions {
    int lloyd_iters_to_add_seed = 5;  // Lloyd iterations before checking oversized clusters
    int lloyd_max_iterations = 100;   // Maximum number of Lloyd iterations
    int random_seed = 42;             // seed for reproducible random picks
};

// ---------------------------------------------------------------------------
// Lloyd namespace: one function per algorithm step.
//
// CC finding lives in ConnectedComponent.h.
// This namespace holds the Lloyd-specific clustering logic, including
// BFS primitives that operate on the CSR graph.
//
// Graph convention everywhere:
//   Gp  -- row pointers   (size G_N + 1)
//   Gi  -- column indices  (size Gp[G_N])
//   G_N -- number of vertices
// ---------------------------------------------------------------------------
namespace Lloyd {

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

// Step 3 -- Initial seed selection (across all CCs).
// For each CC, pick ceil(|CC| / patch_size) random vertices as seeds.
// Output: seeds (populated by the callee).
void initialize_seeds(
    const std::vector<int>& vertex_to_cc,
    const std::vector<int>& num_vertices_per_cc,
    int patch_size,
    int random_seed,
    std::vector<int>& seeds);

// Step 4.3 -- Update cluster centers.
// For each cluster:
//   1. Find boundary vertices (have a neighbour in a different cluster).
//   2. Multi-source BFS inward (restricted to the cluster) from the boundary.
//   3. Pick the new center from the last (deepest) wave.
// Output: seeds is updated in-place with the new centers.
void update_centers(
    const int* Gp, const int* Gi, int G_N,
    const std::vector<int>& node_to_cluster,
    int num_clusters,
    int random_seed,
    std::vector<int>& seeds);

// Step 4.4 -- Split oversized clusters.
// For every cluster whose size exceeds patch_size, add a new seed
// (chosen from the deepest interior wave of that cluster).
// Returns the number of new seeds added (0 means all clusters fit).
void merge_and_split_clusters(
    const int* Gp, const int* Gi, int G_N,
    const std::vector<int>& node_to_cluster,
    const std::vector<int>& cluster_sizes,
    int num_clusters,
    int patch_size,
    int random_seed,
    std::vector<int>& seeds);

} // namespace Lloyd
