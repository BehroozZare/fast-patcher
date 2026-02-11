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
// CC finding lives in ConnectedComponent.h; BFS primitives in BFS.h.
// This namespace holds the Lloyd-specific clustering logic that builds
// on top of those primitives.
//
// Graph convention everywhere:
//   Gp  -- row pointers   (size G_N + 1)
//   Gi  -- column indices  (size Gp[G_N])
//   G_N -- number of vertices
// ---------------------------------------------------------------------------
namespace Lloyd {

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
int split_oversized(
    const int* Gp, const int* Gi, int G_N,
    const std::vector<int>& node_to_cluster,
    const std::vector<int>& cluster_sizes,
    int num_clusters,
    int patch_size,
    int random_seed,
    std::vector<int>& seeds);

} // namespace Lloyd
