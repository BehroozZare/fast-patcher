//
// Created by behrooz on 2025-09-28.
//
#include <igl/read_triangle_mesh.h>
#include <igl/readMESH.h>
#include <igl/cotmatrix.h>
#include <spdlog/spdlog.h>
#include <CLI/CLI.hpp>
#include <Eigen/Core>
#include <chrono>
#include <unsupported/Eigen/SparseExtra>
#include <iostream>
#include <filesystem>
#include <algorithm>
#include <numeric>

#include "remove_diagonal.h"
#include "clusterAPI.h"
#include "polyscope_visualizer.h"

struct CLIArgs
{
    std::string input_mesh;
    std::string output_dir = "./output";
    int patch_size = 512;

    CLIArgs(int argc, char* argv[])
    {
        CLI::App app{"Graph Clustering Demo"};
        app.add_option("-i,--input", input_mesh, "input mesh file (.mesh for tet, .obj/.off/.ply/.stl for triangle)")->required();
        app.add_option("-o,--output", output_dir, "output directory for patch files");
        app.add_option("-z,--patch_size", patch_size, "maximum patch size");

        try {
            app.parse(argc, argv);
        } catch (const CLI::ParseError& e) {
            exit(app.exit(e));
        }
    }
};

// Get file extension in lowercase
std::string get_extension(const std::string& path) {
    std::filesystem::path p(path);
    std::string ext = p.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return ext;
}

// Get mesh name without extension
std::string get_mesh_name(const std::string& path) {
    std::filesystem::path p(path);
    return p.stem().string();
}

// Compute cluster statistics from vertex-to-cluster label vector
struct ClusterStats {
    int num_clusters;
    int min_size;
    int max_size;
    double avg_size;
};

ClusterStats compute_cluster_stats(const std::vector<int>& vertex_to_cluster) {
    ClusterStats stats;

    if (vertex_to_cluster.empty()) {
        stats.num_clusters = 0;
        stats.min_size = 0;
        stats.max_size = 0;
        stats.avg_size = 0.0;
        return stats;
    }

    // Number of clusters = max label + 1
    stats.num_clusters = *std::max_element(vertex_to_cluster.begin(),
                                           vertex_to_cluster.end()) + 1;

    // Count vertices in each cluster
    std::vector<int> sizes(stats.num_clusters, 0);
    for (int label : vertex_to_cluster)
        sizes[label]++;

    stats.min_size = *std::min_element(sizes.begin(), sizes.end());
    stats.max_size = *std::max_element(sizes.begin(), sizes.end());
    int total = static_cast<int>(vertex_to_cluster.size());
    stats.avg_size = static_cast<double>(total) / stats.num_clusters;
    return stats;
}

int main(int argc, char* argv[])
{
    CLIArgs args(argc, argv);

    // Create output directory if it doesn't exist
    std::filesystem::create_directories(args.output_dir);
    
    std::string mesh_name = get_mesh_name(args.input_mesh);
    std::string ext = get_extension(args.input_mesh);
    
    spdlog::info("Loading mesh from: {}", args.input_mesh);
    spdlog::info("Output directory: {}", args.output_dir);
    spdlog::info("Patch size: {}", args.patch_size);

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;  // Elements (triangles or tetrahedra)
    
    // Auto-detect mesh format based on file extension
    bool is_tet_mesh = (ext == ".mesh");
    
    if (is_tet_mesh) {
        // Tet mesh
        Eigen::MatrixXi OF;  // Boundary faces (unused)
        if (!igl::readMESH(args.input_mesh, V, F, OF)) {
            spdlog::error("Failed to read tet mesh: {}", args.input_mesh);
            return 1;
        }
        spdlog::info("Loaded tet mesh: {} vertices, {} tetrahedra", V.rows(), F.rows());
    } else if (ext == ".obj" || ext == ".off" || ext == ".ply" || ext == ".stl") {
        // Triangle mesh
        if (!igl::read_triangle_mesh(args.input_mesh, V, F)) {
            spdlog::error("Failed to read triangle mesh: {}", args.input_mesh);
            return 1;
        }
        spdlog::info("Loaded triangle mesh: {} vertices, {} triangles", V.rows(), F.rows());
    } else {
        spdlog::error("Unsupported mesh format: {}. Supported: .mesh (tet), .obj/.off/.ply/.stl (triangle)", ext);
        return 1;
    }

    int num_vertices = static_cast<int>(V.rows());

    // Compute cotangent Laplacian
    spdlog::info("Computing cotangent matrix...");
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V, F, L);
    spdlog::info("Matrix size: {} x {}, nnz: {}", L.rows(), L.cols(), L.nonZeros());

    // Extract graph without diagonal (undirected adjacency in CSR)
    std::vector<int> Gp, Gi;
    Patcher::remove_diagonal(
        L.rows(), L.outerIndexPtr(), L.innerIndexPtr(), Gp, Gi);
    spdlog::info("Graph: {} vertices, {} edges", num_vertices, Gi.size());

    spdlog::info("========================================");
    spdlog::info("Running Lloyd clustering...");
    spdlog::info("========================================");

    // Time the clustering
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<int> vertex_to_cluster;
    Patcher::create_clusters(
        num_vertices,
        Gp.data(),
        Gi.data(),
        args.patch_size,
        nullptr,
        vertex_to_cluster);
    auto end = std::chrono::high_resolution_clock::now();
    double runtime_ms = std::chrono::duration<double, std::milli>(end - start).count();
    
    // Correctness check: verify size matches number of vertices
    bool size_correct = (static_cast<int>(vertex_to_cluster.size()) == num_vertices);
    
    // Check that all vertices are assigned (no -1 values)
    int unassigned = std::count(vertex_to_cluster.begin(), vertex_to_cluster.end(), -1);
    bool all_assigned = (unassigned == 0);
    
    // Compute statistics
    ClusterStats stats = compute_cluster_stats(vertex_to_cluster);
    
    // Log results
    spdlog::info("  Runtime: {:.2f} ms", runtime_ms);
    spdlog::info("  Clusters: {}, Min: {}, Max: {}, Avg: {:.1f}", 
                    stats.num_clusters, stats.min_size, stats.max_size, stats.avg_size);
    
    if (size_correct && all_assigned) {
        spdlog::info("  Correctness check: PASSED (vertex_to_cluster.size() == {})", num_vertices);
    } else {
        if (!size_correct) {
            spdlog::error("  Correctness check: FAILED (vertex_to_cluster.size() = {}, expected {})", 
                            vertex_to_cluster.size(), num_vertices);
        }
        if (!all_assigned) {
            spdlog::error("  Correctness check: FAILED ({} vertices unassigned)", unassigned);
        }
    }

    spdlog::info("");
    spdlog::info("========================================");
    spdlog::info("Done!");
    spdlog::info("========================================");


    // Visualize the clusters
    PatchVisualizer::visualize_patches(V, F, vertex_to_cluster, is_tet_mesh);

    return 0;
}