#pragma once
#include <vector>
#include <tuple>
#include "bit_array.h"

/**
 * @brief Configuration options for the Lloyd-style clustering algorithm.
 *
 * Controls seed selection strategy, iteration limits, Morton-code
 * quantization resolution, and the random seed for reproducibility.
 */
struct LloydOptions {
    /**
     * @brief Strategy used to choose initial cluster seeds.
     */
    enum class SeedSelectionMethod {
        RANDOM,      ///< Pick seeds uniformly at random within each connected component.
        MORTON_CODE, ///< Pick seeds along a space-filling Morton (Z-order) curve.
        FPS          ///< Pick seeds via farthest point sampling on graph distances.
    };

    SeedSelectionMethod seed_selection_method = SeedSelectionMethod::RANDOM; ///< Seed selection strategy (default: RANDOM).
    int dimension = 3;                       ///< Spatial dimension of vertex positions (used by Morton-code seeding).
    int lloyd_iters_to_add_seed = 1;         ///< Number of Lloyd iterations to run before checking for oversized clusters.
    int lloyd_max_iterations = 100;          ///< Maximum total number of outer Lloyd iterations.
    int morton_code_bits_per_dimension = 10;  ///< Number of quantization bits per axis for Morton-code computation.
    int random_seed = 42;                    ///< RNG seed for reproducible random picks.
};

/**
 * @brief Lloyd-style iterative graph-clustering primitives.
 *
 * This namespace holds the BFS-based building blocks that drive the
 * clustering algorithm: multi-source BFS, seed initialization (random
 * and Morton-code), center updates, and cluster merge/split logic.
 *
 * Connected-component detection lives in ConnectedComponent.h.
 *
 * @par Graph convention (CSR format)
 * | Symbol | Meaning                              |
 * |--------|--------------------------------------|
 * | Gp     | Row pointers (size G_N + 1)          |
 * | Gi     | Column indices (size Gp[G_N])        |
 * | G_N    | Number of vertices                   |
 */
namespace Lloyd {

/**
 * @brief Multi-source BFS on the full CSR graph.
 *
 * Expands level-by-level from every source vertex (those with
 * dist[v] == 0 on entry).  On return every reachable vertex is
 * labelled with the source that first reached it.
 *
 * The BFS is parallelized with OpenMP.
 *
 * @param Gp              Row pointers of the CSR adjacency.
 * @param Gi              Column indices of the CSR adjacency.
 * @param G_N             Number of vertices.
 * @param source_vertices Vector of source vertex indices.
 * @param[in,out] label   On entry, label[s] must be set for each source s.
 *                        On return, label[v] = label of the source that
 *                        first reached v.
 * @param[in,out] dist    On entry, dist[s] = 0 for sources, -1 elsewhere.
 *                        On return, dist[v] = BFS distance from v to its
 *                        nearest source.
 * @param active_mask     Optional bit mask. If non-null, only vertices
 *                        with active_mask->get(v)==true are expanded into.
 *
 * @see restricted_bfs
 */
void multi_source_bfs(
    const int* Gp, const int* Gi, int G_N,
    const std::vector<int>& source_vertices,
    std::vector<int>& label,
    std::vector<int>& dist,
    const BitArray* active_mask = nullptr);

/**
 * @brief Cluster-restricted multi-source BFS.
 *
 * Same traversal logic as multi_source_bfs, but expansion is restricted
 * to vertices that belong to the same cluster as the source.  Returns
 * the vertices reached in the deepest (last) BFS wave.
 *
 * @param Gp                  Row pointers of the CSR adjacency.
 * @param Gi                  Column indices of the CSR adjacency.
 * @param G_N                 Number of vertices.
 * @param node_to_cluster     Cluster assignment for every vertex.
 * @param source_vertices     Boundary or seed vertices to expand from.
 * @param[out] final_wave_vertices  Filled with the vertices of the last
 *                             (deepest) BFS wave.
 * @param[in,out] dist        Distance array (sized G_N, initialised to -1).
 * @param active_mask         Optional bit mask; vertices with bit == 0
 *                            are skipped.
 *
 * @see multi_source_bfs, update_centers
 */
void restricted_bfs(
    const int* Gp, const int* Gi, int G_N,
    const std::vector<int>& node_to_cluster,
    const std::vector<int>& source_vertices,
    std::vector<int>& final_wave_vertices,
    std::vector<int>& dist,
    const BitArray* active_mask = nullptr);

/**
 * @brief Randomly select initial cluster seeds across all connected components.
 *
 * For each connected component, ceil(|CC| / patch_size) vertices are
 * chosen uniformly at random as initial seeds.
 *
 * @param vertex_to_cc         CC label for every vertex.
 * @param num_vertices_per_cc  Number of vertices in each CC.
 * @param patch_size           Target maximum cluster size.
 * @param opt                  Algorithm options (uses random_seed).
 * @param[out] seeds           Populated with the chosen seed vertex indices.
 *
 * @see initialize_seeds_morton_code
 */


void furthest_point_bfs(const int* Gp, const int* Gi, int G_N,
    const int& start_node,
    std::vector<int>& dist,
    int& furthest_distance,
    int& furthest_node);



void initialize_seeds_random(
    const std::vector<int>& vertex_to_cc,
    const std::vector<int>& num_vertices_per_cc,
    int patch_size,
    const LloydOptions& opt,
    std::vector<int>& seeds);

/**
 * @brief Compute the 3D Morton (Z-order) code by bit-interleaving.
 *
 * Interleaves the low-order bits of the three integer coordinates into
 * a single integer whose ordering approximates spatial locality.
 *
 * @param x         Quantized x-coordinate.
 * @param y         Quantized y-coordinate.
 * @param z         Quantized z-coordinate.
 * @param num_bits  Number of low-order bits per coordinate to interleave.
 * @return The interleaved Morton code.
 */
int morton_code(int x, int y, int z, int num_bits);

/**
 * @brief Select initial cluster seeds using a Morton-code space-filling curve.
 *
 * Vertices within each connected component are sorted along the Morton
 * curve (after per-CC bounding-box normalization and quantization).
 * Seeds are then chosen at evenly spaced intervals along the sorted
 * order, yielding spatially well-distributed initial placements.
 *
 * @param vertex_to_cc         CC label for every vertex.
 * @param num_vertices_per_cc  Number of vertices in each CC.
 * @param vertex_positions     3D positions for all vertices (x, y, z).
 * @param patch_size           Target maximum cluster size.
 * @param opt                  Algorithm options (uses morton_code_bits_per_dimension).
 * @param[out] seeds           Populated with the chosen seed vertex indices.
 *
 * @see initialize_seeds_random, morton_code
 */
void initialize_seeds_morton_code(
    const std::vector<int>& vertex_to_cc,
    const std::vector<int>& num_vertices_per_cc,
    const std::vector<std::tuple<double, double, double>>& vertex_positions,
    int patch_size,
    const LloydOptions& opt,
    std::vector<int>& seeds,
    std::vector<int>& sorted_indices);


/**
 * @brief Select initial cluster seeds using a furthest point sampling-based algorithm.
 *
 * Vertices within each connected component are sampled using the furthest point sampling algorithm (FPS).
 *
 * @param Gp                  Row pointers of the CSR adjacency.
 * @param Gi                  Column indices of the CSR adjacency.
 * @param G_N                 Number of vertices.
 * @param vertex_to_cc         CC label for every vertex.
 * @param num_vertices_per_cc  Number of vertices in each CC.
 * @param vertex_positions     3D positions for all vertices (x, y, z).
 * @param patch_size           Target maximum cluster size.
 * @param opt                  Algorithm options (uses random_seed).
 * @param[out] seeds           Populated with the chosen seed vertex indices.
 *
 * @see initialize_seeds_random, initialize_seeds_morton_code
 */
void initialize_seeds_FPS_based(
    const int* Gp, const int* Gi, int G_N,
    const std::vector<int>& vertex_to_cc,
    const std::vector<int>& num_vertices_per_cc,
    const std::vector<std::tuple<double, double, double>>& vertex_positions,
    int patch_size,
    const LloydOptions& opt,
    std::vector<int>& seeds);

/**
 * @brief Update cluster centres to the deepest interior vertex.
 *
 * For each cluster:
 *   1. Identify boundary vertices (those with at least one neighbour in
 *      a different cluster).
 *   2. Run a restricted BFS inward from the boundary (within the cluster).
 *   3. Select a new seed from the last (deepest) BFS wave.
 *
 * Seeds are updated in-place.
 *
 * @param Gp              Row pointers of the CSR adjacency.
 * @param Gi              Column indices of the CSR adjacency.
 * @param G_N             Number of vertices.
 * @param node_to_cluster Current cluster assignment for every vertex.
 * @param num_clusters    Total number of clusters.
 * @param opt             Algorithm options (uses random_seed for tie-breaking).
 * @param[in,out] seeds   On entry, current seed vertices; on return,
 *                        updated to the new centres.
 * @param active_mask     Optional bit mask; clusters whose seed has
 *                        bit == 0 are skipped.
 *
 * @see restricted_bfs
 */
void update_centers(
    const int* Gp, const int* Gi, int G_N,
    const std::vector<int>& node_to_cluster,
    int num_clusters,
    const LloydOptions& opt,
    std::vector<int>& seeds,
    const BitArray* active_mask = nullptr);

/**
 * @brief Merge undersized clusters and split oversized ones.
 *
 * **Merge phase**: Adjacent cluster pairs whose combined size does not
 * exceed @p patch_size are greedily merged, prioritising the pair with
 * the smallest combined size.  Merging removes one seed (the absorbed
 * cluster) so the next BFS round reassigns its vertices.
 *
 * **Split phase**: For every remaining cluster whose size exceeds
 * @p patch_size, a new seed is injected at a neighbour of the existing
 * seed that belongs to the same cluster.
 *
 * @param Gp              Row pointers of the CSR adjacency.
 * @param Gi              Column indices of the CSR adjacency.
 * @param G_N             Number of vertices.
 * @param node_to_cluster Current cluster assignment for every vertex.
 * @param cluster_sizes   Number of vertices in each cluster.
 * @param num_clusters    Total number of clusters before this call.
 * @param patch_size      Maximum allowed cluster size.
 * @param opt             Algorithm options.
 * @param[in,out] seeds   Seed list; absorbed seeds are removed, new
 *                        seeds are appended for split clusters.
 * @param[out] cluster_was_merged  Resized to @p num_clusters.  Entry
 *                        [i] is true iff cluster i was absorbed into a
 *                        neighbour during the merge phase.
 */
void merge_and_split_clusters(
    const int* Gp, const int* Gi, int G_N,
    std::vector<int>& node_to_cluster,
    const std::vector<int>& cluster_sizes,
    int num_clusters,
    int patch_size,
    const LloydOptions& opt,
    std::vector<int>& seeds,
    std::vector<bool>& cluster_was_merged);

} // namespace Lloyd
