#pragma once

#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <array>

#include <Eigen/Core>
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/volume_mesh.h"

namespace PatchVisualizer {

// Convert HSL to RGB (h in [0,360), s,l in [0,1]) → RGB in [0,1]
inline std::array<double, 3> hsl_to_rgb(double h, double s, double l) {
    double c = (1.0 - std::abs(2.0 * l - 1.0)) * s;
    double hp = h / 60.0;
    double x = c * (1.0 - std::abs(std::fmod(hp, 2.0) - 1.0));
    double r1 = 0, g1 = 0, b1 = 0;

    if      (hp < 1) { r1 = c; g1 = x; }
    else if (hp < 2) { r1 = x; g1 = c; }
    else if (hp < 3) { g1 = c; b1 = x; }
    else if (hp < 4) { g1 = x; b1 = c; }
    else if (hp < 5) { r1 = x; b1 = c; }
    else             { r1 = c; b1 = x; }

    double m = l - c / 2.0;
    return {r1 + m, g1 + m, b1 + m};
}

// Generate `n` visually distinct colors using evenly spaced hues with jitter.
inline Eigen::MatrixXd generate_cluster_colors(int n, int seed = 42) {
    Eigen::MatrixXd colors(n, 3);
    if (n == 0) return colors;
    if (n == 1) {
        colors.row(0) << 0.267, 0.467, 0.667; // pleasant blue
        return colors;
    }

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> jitter(-15.0, 15.0);
    std::uniform_real_distribution<double> sat_dist(0.55, 0.85);
    std::uniform_real_distribution<double> lit_dist(0.45, 0.65);

    double hue_step = 360.0 / n;
    for (int i = 0; i < n; ++i) {
        double hue = std::fmod(i * hue_step + jitter(rng) + 360.0, 360.0);
        double sat = sat_dist(rng);
        double lit = lit_dist(rng);
        auto rgb = hsl_to_rgb(hue, sat, lit);
        colors.row(i) << rgb[0], rgb[1], rgb[2];
    }
    return colors;
}

// Visualize mesh patches using Polyscope.
//
//   V                - Nx3 vertex positions
//   F                - element connectivity (Fx3 triangles or Fx4 tetrahedra)
//   vertex_to_cluster - per-vertex cluster label (size N)
//   is_tet_mesh      - true for volumetric (tet) meshes
//
inline void visualize_patches(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const std::vector<int>& vertex_to_cluster,
    bool is_tet_mesh)
{
    // --- Polyscope init ---
    polyscope::options::programName = "Fast Patcher – Patch Viewer";
    polyscope::init();

    int num_clusters = *std::max_element(vertex_to_cluster.begin(),
                                         vertex_to_cluster.end()) + 1;

    // Build per-vertex color matrix from cluster labels
    Eigen::MatrixXd cluster_colors = generate_cluster_colors(num_clusters);
    int num_verts = static_cast<int>(V.rows());
    Eigen::MatrixXd vertex_colors(num_verts, 3);
    for (int i = 0; i < num_verts; ++i) {
        vertex_colors.row(i) = cluster_colors.row(vertex_to_cluster[i]);
    }

    if (is_tet_mesh) {
        // Register as volume mesh (tet)
        auto* vm = polyscope::registerVolumeMesh("mesh", V, F);
        vm->addVertexColorQuantity("patch", vertex_colors)->setEnabled(true);
    } else {
        // Register as surface mesh (triangle)
        auto* sm = polyscope::registerSurfaceMesh("mesh", V, F);
        sm->addVertexColorQuantity("patch", vertex_colors)->setEnabled(true);
    }

    polyscope::show();
}

} // namespace PatchVisualizer
