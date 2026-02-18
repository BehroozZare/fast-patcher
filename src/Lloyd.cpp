#include "Lloyd.h"
#include "bit_array.h"
#include <cmath>
#include <random>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <limits>
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
    std::vector<int>& dist,
    const BitArray* active_mask)
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
                    if (active_mask && !active_mask->get(u)) continue;
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
    std::vector<int>& dist,
    const BitArray* active_mask)
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
                if (active_mask && !active_mask->get(u)) continue;
                if (node_to_cluster[u] != cluster_id) continue;
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


void furthest_point_bfs(const int* Gp, const int* Gi, int G_N,
    const int& start_node,
    std::vector<int>& dist,
    int& furthest_distance,
    int& furthest_node){

    Lloyd::BitArray visited(G_N);
    visited.set(start_node);
    dist[start_node] = 0;
    std::vector<int> current_queue{start_node};
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
                        if(dist[u] > dist[v] + 1) {
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

    // Post-BFS: scan all vertices to find the true argmax of min-distance
    furthest_distance = 0;
    furthest_node = start_node;
    for (int v = 0; v < G_N; v++) {
        if (dist[v] != std::numeric_limits<int>::max() && dist[v] > furthest_distance) {
            furthest_distance = dist[v];
            furthest_node = v;
        }
    }
}

// -----------------------------------------------------------------------
// Step 3 -- Initial seed selection (across all CCs).
//
// For each CC, pick ceil(|CC| / patch_size) random vertices as seeds.
// -----------------------------------------------------------------------
void initialize_seeds_random(
    const std::vector<int>& vertex_to_cc,
    const std::vector<int>& num_vertices_per_cc,
    int patch_size,
    const LloydOptions& opt,
    std::vector<int>& seeds)
{
    seeds.clear();
    //Set the random seed
    std::mt19937 rng(opt.random_seed);
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
    const LloydOptions& opt,
    std::vector<int>& seeds,
    const BitArray* active_mask)
{
    //Find the center per each cluster
    //Step 0: Find all the boundary vertices
    std::vector<std::vector<int>> per_cluster_boundary_vertices(num_clusters);
    for(int v = 0; v < G_N; v++) {
        if (active_mask && !active_mask->get(v)) continue;
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
    #pragma omp parallel for
    for(int i = 0; i < num_clusters; i++) {
        if (active_mask && !active_mask->get(seeds[i])) continue;
        std::vector<int> final_wave_vertices;
        Lloyd::restricted_bfs(Gp, Gi, G_N, node_to_cluster, per_cluster_boundary_vertices[i], final_wave_vertices, dist, active_mask);
        //Pick a random new seed from the final wave vertices
        std::mt19937 rng(opt.random_seed + i);
        std::uniform_int_distribution<int> dist_rng(0, final_wave_vertices.size() - 1);
        seeds[i] = final_wave_vertices[dist_rng(rng)];
    }
}


int morton_code(int x, int y, int z, int morton_code_bits_per_dimension) {
    int morton_code = 0;
    for(int i = 0; i < morton_code_bits_per_dimension; i++) {
        morton_code |= (x & (1 << i)) << (2 * i);
        morton_code |= (y & (1 << i)) << (2 * i + 1);
        morton_code |= (z & (1 << i)) << (2 * i + 2);
    }
    return morton_code;
}

// -----------------------------------------------------------------------
// Step 3 -- Initial seed selection (across all CCs).
//
// For each CC, pick ceil(|CC| / patch_size) random vertices as seeds.
// -----------------------------------------------------------------------
void initialize_seeds_morton_code(
    const std::vector<int>& vertex_to_cc,
    const std::vector<int>& num_vertices_per_cc,
    const std::vector<std::tuple<double, double, double>>& vertex_positions,
    int patch_size,
    const LloydOptions& opt,
    std::vector<int>& seeds,
    std::vector<int>& sorted_indices)
{
    seeds.clear();
    int morton_range = (1 << opt.morton_code_bits_per_dimension) - 1;
    // Group vertices by their CC
    std::vector<std::vector<int>> cc_vertices(num_vertices_per_cc.size());
    for (int i = 0; i < static_cast<int>(vertex_to_cc.size()); i++) {
        cc_vertices[vertex_to_cc[i]].push_back(i);
    }

    // Map seed indices to actual vertex ids
    for (int cc = 0; cc < static_cast<int>(cc_vertices.size()); cc++) {
        int num_seeds_for_cc = std::ceil(static_cast<double>(num_vertices_per_cc[cc]) / static_cast<double>(patch_size));

        // Compute per-CC bounding box
        double cc_min_x = std::get<0>(vertex_positions[cc_vertices[cc][0]]);
        double cc_min_y = std::get<1>(vertex_positions[cc_vertices[cc][0]]);
        double cc_min_z = std::get<2>(vertex_positions[cc_vertices[cc][0]]);
        double cc_max_x = cc_min_x, cc_max_y = cc_min_y, cc_max_z = cc_min_z;
        for (int i = 1; i < num_vertices_per_cc[cc]; i++) {
            int vid = cc_vertices[cc][i];
            double x = std::get<0>(vertex_positions[vid]);
            double y = std::get<1>(vertex_positions[vid]);
            double z = std::get<2>(vertex_positions[vid]);
            cc_min_x = std::min(cc_min_x, x); cc_max_x = std::max(cc_max_x, x);
            cc_min_y = std::min(cc_min_y, y); cc_max_y = std::max(cc_max_y, y);
            cc_min_z = std::min(cc_min_z, z); cc_max_z = std::max(cc_max_z, z);
        }
        double ext_x = cc_max_x - cc_min_x; if (ext_x == 0.0) ext_x = 1.0;
        double ext_y = cc_max_y - cc_min_y; if (ext_y == 0.0) ext_y = 1.0;
        double ext_z = cc_max_z - cc_min_z; if (ext_z == 0.0) ext_z = 1.0;

        // Quantize positions and compute morton codes
        std::vector<int> morton_codes(num_vertices_per_cc[cc]);
        for (int i = 0; i < num_vertices_per_cc[cc]; i++) {
            int vertex_id = cc_vertices[cc][i];
            int qx = static_cast<int>((std::get<0>(vertex_positions[vertex_id]) - cc_min_x) / ext_x * morton_range);
            int qy = static_cast<int>((std::get<1>(vertex_positions[vertex_id]) - cc_min_y) / ext_y * morton_range);
            int qz = static_cast<int>((std::get<2>(vertex_positions[vertex_id]) - cc_min_z) / ext_z * morton_range);
            morton_codes[i] = morton_code(qx, qy, qz, opt.morton_code_bits_per_dimension);
        }
        // Sort vertex indices by their morton code (preserves mapping to cc_vertices)
        sorted_indices.resize(num_vertices_per_cc[cc]);
        std::iota(sorted_indices.begin(), sorted_indices.end(), 0);
        std::sort(sorted_indices.begin(), sorted_indices.end(), [&morton_codes](int a, int b) {
            return morton_codes[a] < morton_codes[b];
        });

        if (num_vertices_per_cc[cc] <= num_seeds_for_cc) {
            // Fewer vertices than seeds: pick every vertex as a seed
            for (int i = 0; i < num_vertices_per_cc[cc]; i++) {
                seeds.push_back(cc_vertices[cc][sorted_indices[i]]);
            }
        } else {
            // Divide sorted order into num_seeds_for_cc equal parts, pick the middle of each
            int part_size = num_vertices_per_cc[cc] / num_seeds_for_cc;
            for (int i = 0; i < num_seeds_for_cc; i++) {
                int mid = i * part_size + part_size / 2;
                seeds.push_back(cc_vertices[cc][sorted_indices[mid]]);
            }
        }
    }

}

void initialize_seeds_FPS_based(
    const int* Gp, const int* Gi, int G_N,
    const std::vector<int>& vertex_to_cc,
    const std::vector<int>& num_vertices_per_cc,
    const std::vector<std::tuple<double, double, double>>& vertex_positions,
    int patch_size,
    const LloydOptions& opt,
    std::vector<int>& seeds)
{
    seeds.clear();
    for(int i = 0; i < num_vertices_per_cc.size(); i++) {
        //Pick seeds per each CC using the FPS algorithm
        //Step 0: Find a single node from the CC to start the FPS algorithm
        int start_node = -1;
        for(int j = 0; j < G_N; j++) {
            if(vertex_to_cc[j] == i) {
                start_node = j;
                break;
            }
        }
        //Step 1: Find then umber of seeds to pick for this CC
        int num_seeds_for_cc = std::ceil(static_cast<double>(num_vertices_per_cc[i]) / static_cast<double>(patch_size));
        //Step 2: Iteratively find the seeds by finding the farthest node from the current seed
        int furthest_distance = 0;
        int furthest_node = start_node;
        std::vector<int> dist(G_N, std::numeric_limits<int>::max());
        for(int j = 0; j < num_seeds_for_cc; j++) {
            Lloyd::furthest_point_bfs(Gp, Gi, G_N, furthest_node, dist, furthest_distance, furthest_node);
            seeds.push_back(furthest_node);
        }
    }
}
// -----------------------------------------------------------------------
// Step 4.4 -- Merge undersized clusters, then split oversized ones.
// -----------------------------------------------------------------------
void merge_and_split_clusters(
    const int* Gp, const int* Gi, int G_N,
    std::vector<int>& node_to_cluster,
    const std::vector<int>& cluster_sizes,
    int num_clusters,
    int patch_size,
    const LloydOptions& opt,
    std::vector<int>& seeds,
    std::vector<bool>& cluster_was_merged)
{
    // ── Merge phase ─────────────────────────────────────────────────
    // 1. Build unique cluster-pair adjacency by scanning all edges.
    std::vector<std::pair<int,int>> adj_pairs;
    for (int v = 0; v < G_N; v++) {
        int cv = node_to_cluster[v];
        for (int j = Gp[v]; j < Gp[v + 1]; j++) {
            int cu = node_to_cluster[Gi[j]];
            if (cv < cu)
                adj_pairs.push_back({cv, cu});
        }
    }
    std::sort(adj_pairs.begin(), adj_pairs.end());
    adj_pairs.erase(std::unique(adj_pairs.begin(), adj_pairs.end()),
                    adj_pairs.end());

    // 2. Collect merge candidates (combined size <= patch_size).
    struct MergeCandidate { int combined; int a; int b; };
    std::vector<MergeCandidate> candidates;
    for (auto& [a, b] : adj_pairs) {
        int combined = cluster_sizes[a] + cluster_sizes[b];
        if (combined <= patch_size)
            candidates.push_back({combined, a, b});
    }

    // 3. Sort by combined size ascending (smallest merges first).
    std::sort(candidates.begin(), candidates.end(),
        [](const MergeCandidate& x, const MergeCandidate& y) {
            return x.combined < y.combined;
        });

    // 4. Greedy matching -- each cluster participates in at most one merge.
    cluster_was_merged.assign(num_clusters, false);
    std::vector<bool> involved(num_clusters, false);
    std::vector<int> cluster_mapped_to(num_clusters, -1);
    int merge_count = 0;
    for (auto& c : candidates) {
        if (!involved[c.a] && !involved[c.b]) {
            involved[c.a] = true;
            involved[c.b] = true;
            cluster_was_merged[c.b] = true;
            merge_count++;
            cluster_mapped_to[c.b] = c.a;
        }
    }

    // 5. Remove absorbed seeds via swap-and-pop (descending order).
    if (merge_count > 0) {
        std::vector<int> to_remove;
        to_remove.reserve(merge_count);
        for (int i = 0; i < num_clusters; i++)
            if (cluster_was_merged[i])
                to_remove.push_back(i);
        std::sort(to_remove.rbegin(), to_remove.rend());

        std::vector<int> slot_occupant(num_clusters);
        std::iota(slot_occupant.begin(), slot_occupant.end(), 0);

        for (int idx : to_remove) {
            int last = static_cast<int>(seeds.size()) - 1;
            if (idx != last) {
                seeds[idx] = seeds[last];
                int occupant = slot_occupant[last];
                slot_occupant[idx] = occupant;
                cluster_mapped_to[occupant] = idx;
            }
            seeds.pop_back();
        }
    }

    #pragma omp parallel for
    for (int v = 0; v < G_N; v++) {
        int cluster_id = node_to_cluster[v];
        if (cluster_mapped_to[cluster_id] != -1) {
            node_to_cluster[v] = cluster_mapped_to[cluster_id];
        }
    }

    spdlog::info("  merge_and_split: merged {} cluster pairs", merge_count);

    // ── Split phase ─────────────────────────────────────────────────
    // Add a new seed at a same-cluster neighbour of each oversized
    // cluster's seed.  Only consider seeds that survived the merge.
    std::vector<int> added_seeds;
    for (int seed : seeds) {
        int cluster_id = node_to_cluster[seed];
        if (cluster_sizes[cluster_id] < patch_size)
            continue;
        for (int j = Gp[seed]; j < Gp[seed + 1]; j++) {
            int neighbor = Gi[j];
            if (node_to_cluster[neighbor] == cluster_id) {
                added_seeds.push_back(neighbor);
                break;
            }
        }
    }
    seeds.insert(seeds.end(), added_seeds.begin(), added_seeds.end());

    spdlog::info("  merge_and_split: split {} oversized clusters", added_seeds.size());
}

} // namespace Lloyd
