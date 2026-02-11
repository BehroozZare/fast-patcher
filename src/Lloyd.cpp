#include "Lloyd.h"
// #include "BFS.h"  -- uncomment when implementing the TODOs below
#include <cassert>
// When implementing the TODOs below you will likely need:
//   #include <algorithm>  -- for std::shuffle, std::min_element, etc.
//   #include <random>     -- for std::mt19937

namespace Lloyd {

// -----------------------------------------------------------------------
// Step 3 -- Initial seed selection (across all CCs).
//
// For each CC, pick ceil(|CC| / patch_size) random vertices as seeds.
// -----------------------------------------------------------------------
void initialize_seeds(
    const std::vector<int>& vertex_to_cc,
    const std::vector<int>& num_vertices_per_cc,
    int patch_size,
    int random_seed,
    std::vector<int>& seeds)
{
    seeds.clear();

    //Generate random numbers based on the number of vertices in each CC
    std::vector<std::vector<int>> per_cc_seeds_ids(num_vertices_per_cc.size());
    for(int i = 0; i < num_vertices_per_cc.size(); i++) {
        //First: Compute how many seeds we need for this cc
        int num_seeds_for_cc = std::ceil(static_cast<double>(num_vertices_per_cc[i]) / static_cast<double>(patch_size));
        int num_vertices_in_cc = num_vertices_per_cc[i];
        for(int j = 0; j < num_seeds_for_cc; j++) {
            int random_seed = std::rand() % num_vertices_in_cc;
            per_cc_seeds_ids[i].push_back(random_seed);
        }
    }
    //Now we have the seeds for each CC, we need to first extract the seed ids from the vertex_to_cc vector
    //First create a pair of vertex, cc_id
    std::vector<std::pair<int, int>> vertex_to_cc_pair;
    for(int i = 0; i < vertex_to_cc.size(); i++) {
        vertex_to_cc_pair.push_back(std::make_pair(i, vertex_to_cc[i]));
    }
    //Now sort the vector based on the cc_id
    std::sort(vertex_to_cc_pair.begin(), vertex_to_cc_pair.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
        return a.second < b.second;
    });
    
    std::vector<int> cumulative_num_vertices_per_cc(num_vertices_per_cc.size() + 1, 0);
    for(int i = 0; i < vertex_to_cc_pair.size(); i++) {
        cumulative_num_vertices_per_cc[i + 1] = cumulative_num_vertices_per_cc[i] + num_vertices_per_cc[i];
    }

    //Now we can extract the seeds for each CC
    for(int cc = 0; cc < num_vertices_per_cc.size(); cc++) {
        for(int seed_id : per_cc_seeds_ids[cc]) {
            int offset = cumulative_num_vertices_per_cc[cc];
            int seed_vertex = seed_id + offset;
            seeds.push_back(vertex_to_cc_pair[seed_vertex].first);
        }
    }
}

// -----------------------------------------------------------------------
// Step 4.3 -- Update cluster centers.
//
// For each cluster C:
//   1. Identify boundary vertices of C (vertices with >= 1 neighbour
//      outside C).
//   2. Multi-source BFS from the boundary, restricted to vertices in C.
//   3. The last wave reached = the deepest interior vertices.
//   4. Pick the new center from that last wave (random or min-id).
// -----------------------------------------------------------------------
void update_centers(
    const int* Gp, const int* Gi, int G_N,
    const std::vector<int>& node_to_cluster,
    int num_clusters,
    int random_seed,
    std::vector<int>& seeds)
{
    // TODO: implement center update using BFS::restricted_bfs.
    //  - For each cluster i in [0, num_clusters):
    //      1. Collect boundary vertices: v in cluster i where at least one
    //         neighbour (via Gp/Gi) has node_to_cluster[neighbour] != i.
    //      2. Prepare label/dist arrays (sized G_N, init to -1).
    //         Set label[b] = 0, dist[b] = 0 for each boundary vertex b.
    //      3. Call BFS::restricted_bfs(Gp, Gi, G_N,
    //                                  node_to_cluster, i, label, dist).
    //      4. Find the maximum distance in dist[] (the deepest wave).
    //         Collect W_last = vertices with that max distance.
    //      5. Pick seeds[i] = a vertex from W_last.
    //         (Use random_seed + i to seed the random choice, or just
    //          pick the vertex with the smallest id for determinism.)
}

// -----------------------------------------------------------------------
// Step 4.4 -- Split oversized clusters.
//
// After some Lloyd iterations, any cluster larger than patch_size gets
// an extra seed injected at its deepest interior vertex.
// -----------------------------------------------------------------------
int split_oversized(
    const int* Gp, const int* Gi, int G_N,
    const std::vector<int>& node_to_cluster,
    const std::vector<int>& cluster_sizes,
    int num_clusters,
    int patch_size,
    int random_seed,
    std::vector<int>& seeds)
{
    int added = 0;

    // TODO: implement oversized-cluster splitting using BFS::restricted_bfs.
    //  - cluster_sizes[i] is already provided (number of vertices in cluster i).
    //  - For every cluster i with cluster_sizes[i] > patch_size:
    //      1. Find its boundary vertices (same as in update_centers).
    //      2. Prepare label/dist arrays (sized G_N, init to -1).
    //         Set label[b] = 0, dist[b] = 0 for each boundary vertex b.
    //      3. Call BFS::restricted_bfs(Gp, Gi, G_N,
    //                                  node_to_cluster, i, label, dist).
    //      4. Find W_last (vertices with max distance).
    //         Pick a new seed from W_last, append it to seeds[].
    //      5. Increment added.
    //  - Return added (0 means convergence: all clusters fit).

    return added;
}

} // namespace Lloyd
