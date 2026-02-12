#include "Lloyd.h"
#include "bit_array.h"
#include <cmath>
#include <random>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <omp.h>
#include <spdlog/spdlog.h>

namespace Lloyd {

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
    Lloyd::BitArray visited(G_N);
    for (int s : source_vertices)
        visited.set(s);

    std::vector<int> current_queue(source_vertices.begin(), source_vertices.end());
    std::vector<int> next_queue;

    #pragma omp parallel
    {
        std::vector<int> local_next_queue;
        while (!current_queue.empty()) {
            #pragma omp for nowait
            for (int i = 0; i < (int)current_queue.size(); i++) {
                int v = current_queue[i];
                for (int j = Gp[v]; j < Gp[v + 1]; j++) {
                    int u = Gi[j];
                    if (!visited.get(u)) {
                        visited.set(u);
                        if(dist[u] == -1) {
                            label[u] = label[v];
                            dist[u]  = dist[v] + 1;
                            local_next_queue.push_back(u);
                        }
                    }
                }
            }
            #pragma omp critical
            {
                next_queue.insert(next_queue.end(),
                    local_next_queue.begin(), local_next_queue.end());
                local_next_queue.clear();
            }
            #pragma omp barrier
            #pragma omp single
            {
                std::swap(current_queue, next_queue);
                next_queue.clear();
            }
        }
    }
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
    const std::vector<int>& node_to_cluster,
    const std::vector<int>& source_vertices,
    std::vector<int>& final_wave_vertices,
    std::vector<int>& dist)
{

    //Initialise visited set from the sources.
    Lloyd::BitArray visited(G_N);
    for (int s : source_vertices) {
        visited.set(s);
        dist[s] = 0;
    }

    //BFS level by level (same OpenMP structure as multi_source_bfs).
    std::vector<int> current_queue(source_vertices.begin(), source_vertices.end());
    std::vector<int> next_queue;
    while (!current_queue.empty()) {
        for (int i = 0; i < (int)current_queue.size(); i++) {
            int v = current_queue[i];
            int cluster_id = node_to_cluster[v];
            for (int j = Gp[v]; j < Gp[v + 1]; j++) {
                int u = Gi[j];
                if (node_to_cluster[u] != cluster_id) continue; //Skip if the neighbor is not in the same cluster
                if (!visited.get(u)) {
                    visited.set(u);
                    if (dist[u] == -1) {
                        dist[u]  = dist[v] + 1;
                        next_queue.push_back(u);
                    }
                }
            }
        }
        if(next_queue.empty()) {
            break;
        }
        std::swap(current_queue, next_queue);
        next_queue.clear();
    }
    std::swap(current_queue, final_wave_vertices);
    current_queue.clear();
}

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
    //Set the random seed
    std::mt19937 rng(random_seed);
    //Generate random numbers based on the number of vertices in each CC
    std::vector<std::vector<int>> per_cc_seeds_ids(num_vertices_per_cc.size());
    for(int i = 0; i < num_vertices_per_cc.size(); i++) {
        //First: Compute how many seeds we need for this cc
        int num_seeds_for_cc = std::ceil(static_cast<double>(num_vertices_per_cc[i]) / static_cast<double>(patch_size));
        int num_vertices_in_cc = num_vertices_per_cc[i];

        // Build index list [0, 1, ..., num_vertices_in_cc - 1]
        std::vector<int> indices(num_vertices_in_cc);
        std::iota(indices.begin(), indices.end(), 0);

        // Shuffle and pick the first num_seeds_for_cc unique indices
        std::shuffle(indices.begin(), indices.end(), rng);
        per_cc_seeds_ids[i].assign(indices.begin(), indices.begin() + num_seeds_for_cc);
    }
    // Group vertices by their CC
    std::vector<std::vector<int>> cc_vertices(num_vertices_per_cc.size());
    for (int i = 0; i < static_cast<int>(vertex_to_cc.size()); i++) {
        cc_vertices[vertex_to_cc[i]].push_back(i);
    }

    // Map seed indices to actual vertex ids
    for (int cc = 0; cc < static_cast<int>(num_vertices_per_cc.size()); cc++) {
        for (int seed_id : per_cc_seeds_ids[cc]) {
            seeds.push_back(cc_vertices[cc][seed_id]);
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
    //Find the center per each cluster
    //Step 0: Find all the boundary vertices
    std::vector<std::vector<int>> per_cluster_boundary_vertices(num_clusters);
    for(int v = 0; v < G_N; v++) {
        for(int j = Gp[v]; j < Gp[v + 1]; j++) {
            int nbr_v = Gi[j];
            if(node_to_cluster[nbr_v] != node_to_cluster[v]) {
                per_cluster_boundary_vertices[node_to_cluster[v]].push_back(v);
                break;
            }
        }
    }

    //Step 1: Run the restricted BFS from the boundary vertices
    //only nbrs with the same cluster id are considered for each boundary vertex
    std::vector<int> dist(G_N, -1);
    seeds.clear();
    for(int i = 0; i < num_clusters; i++) {
        std::vector<int> final_wave_vertices;
        Lloyd::restricted_bfs(Gp, Gi, G_N, node_to_cluster, per_cluster_boundary_vertices[i], final_wave_vertices, dist);
        //Pick a random new seed from the final wave vertices
        std::mt19937 rng(random_seed + i);
        std::uniform_int_distribution<int> dist_rng(0, final_wave_vertices.size() - 1);
        int new_seed = final_wave_vertices[dist_rng(rng)];
        seeds.push_back(new_seed);
    }

#ifndef NDEBUG
    for (int i = 0; i < G_N; i++) {
        assert(dist[i] != -1);
    }
#endif
}

// -----------------------------------------------------------------------
// Step 4.4 -- Split oversized clusters.
//
// After some Lloyd iterations, any cluster larger than patch_size gets
// an extra seed injected at its deepest interior vertex.
// -----------------------------------------------------------------------
void merge_and_split_clusters(
    const int* Gp, const int* Gi, int G_N,
    const std::vector<int>& node_to_cluster,
    const std::vector<int>& cluster_sizes,
    int num_clusters,
    int patch_size,
    int random_seed,
    std::vector<int>& seeds)
{
    int split_clusters = 0;
    int merge_clusters = 0;


    //Go through the seeds, and add a seed to neighbors of that seed
    //if the neighbor is in the same cluster and 
    //the cluster size is more than patch_size
    //then add a new seed to the neighbor
    std::vector<int> added_seeds;
    for(int i = 0; i < seeds.size(); i++) {
        int seed = seeds[i];
        int cluster_id = node_to_cluster[seed];
        if (cluster_sizes[cluster_id] < patch_size) {
            continue;
        }
        for(int j = Gp[seed]; j < Gp[seed + 1]; j++) {
            int neighbor = Gi[j];
            if(node_to_cluster[neighbor] == cluster_id) {
                added_seeds.push_back(neighbor);
                break;
            }
        }
    }
    //Merge added seeds into the seeds
    seeds.insert(seeds.end(), added_seeds.begin(), added_seeds.end());
}

} // namespace Lloyd
