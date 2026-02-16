# fast-patcher

Fast Lloyd-style graph clustering for mesh patching.

## Overview

**fast-patcher** partitions an undirected graph (stored in Compressed Sparse
Row format) into bounded-size clusters using an iterative Lloyd refinement
scheme.  It is designed for mesh patching workflows where a triangle or
tetrahedral mesh must be decomposed into contiguous patches of at most
`patch_size` vertices.

Key features:

- **Bounded cluster size** -- every cluster is guaranteed to contain at most
  `patch_size` vertices after convergence.
- **Two seeding strategies** -- uniform random or spatially-aware Morton
  (Z-order) curve seeding for better initial placement.
- **OpenMP parallelism** -- the multi-source BFS kernel is parallelized for
  multi-core machines.
- **Header-only API** -- a single call to `Patcher::create_clusters()` runs
  the entire pipeline.

## Algorithm

The clustering pipeline executed by `Patcher::create_clusters()` consists of
the following steps:

1. **Connected-component detection** -- label every vertex with its CC id so
   that clustering never spans disconnected subgraphs.
2. **Seed initialization** -- for each CC, select
   `ceil(|CC| / patch_size)` initial seeds using one of:
   - *Random* -- vertices chosen uniformly at random.
   - *Morton-code* -- vertices sorted along a 3D Z-order curve and sampled
     at equal intervals for spatially distributed placement.
3. **Multi-source BFS assignment** -- expand level-by-level from all seeds
   simultaneously; each vertex is assigned to the cluster of the seed that
   reaches it first.
4. **Centre update** -- for each cluster, identify boundary vertices, run an
   inward restricted BFS, and move the seed to the deepest interior vertex.
5. **Merge / split** -- if any cluster exceeds `patch_size`, inject a new
   seed inside it and mark the cluster active for the next round.
6. **Repeat** steps 3--5 until every cluster fits within the size bound or
   the maximum iteration count is reached.

## API Usage

```cpp
#include "clusterAPI.h"

// CSR adjacency of your graph
int n = ...;
std::vector<int> Gp = ...;  // row pointers  (size n+1)
std::vector<int> Gi = ...;  // column indices (size Gp[n])

// Optional: 3D vertex positions for Morton-code seeding
std::vector<std::tuple<double,double,double>> positions = ...;

LloydOptions opts;
opts.seed_selection_method = LloydOptions::SeedSelectionMethod::MORTON_CODE;
opts.lloyd_iters_to_add_seed = 2;

std::vector<int> vertex_to_cluster;
Patcher::create_clusters(
    n, Gp.data(), Gi.data(),
    /*patch_size=*/512,
    &opts,
    vertex_to_cluster,
    &positions);

// vertex_to_cluster[v] is now the cluster id for vertex v
```

Pass `nullptr` for `opt` to use default settings, and `nullptr` for
`vertex_positions` to skip Morton-code seeding (falls back to random).

## Building

### Prerequisites

- C++17 compiler with OpenMP support
- CMake >= 3.16

All other dependencies are fetched automatically via CMake `FetchContent`
recipes located in `cmake/recipes/`:

| Dependency | Purpose                          |
|------------|----------------------------------|
| Eigen      | Matrix types for mesh I/O        |
| libigl     | Mesh reading and Laplacian       |
| spdlog     | Logging                          |
| CLI11      | Command-line argument parsing    |
| Polyscope  | Interactive 3D visualisation     |
| OpenMP     | Parallel BFS                     |

### Build Steps

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

## Demo

The repository ships a demo executable `test_cluster` that reads a mesh,
builds a CSR graph from the cotangent Laplacian, runs the clustering, and
visualises the result with Polyscope.

```bash
./build/test_cluster -i mesh.obj -z 512 -s morton -l 2
```

### CLI Flags

| Flag | Long form        | Description                                         | Default     |
|------|------------------|-----------------------------------------------------|-------------|
| `-i` | `--input`       | Input mesh file (`.mesh` for tet, `.obj`/`.off`/`.ply`/`.stl` for triangle) | *required* |
| `-o` | `--output`      | Output directory for patch files                    | `./output`  |
| `-z` | `--patch_size`  | Maximum number of vertices per cluster              | `512`       |
| `-s` | `--seed`        | Seed selection method: `random` or `morton`          | `morton`    |
| `-l` | `--lloyd_iters` | Lloyd iterations before checking oversized clusters  | `2`         |

## Configuration

All tuneable parameters live in the `LloydOptions` struct:

| Field                          | Type                  | Default    | Description                                              |
|--------------------------------|-----------------------|------------|----------------------------------------------------------|
| `seed_selection_method`        | `SeedSelectionMethod` | `RANDOM`   | Strategy for choosing initial seeds                      |
| `dimension`                    | `int`                 | `3`        | Spatial dimension (used by Morton-code seeding)           |
| `lloyd_iters_to_add_seed`      | `int`                 | `1`        | Lloyd iterations between oversized-cluster checks         |
| `lloyd_max_iterations`         | `int`                 | `100`      | Maximum total outer iterations                           |
| `morton_code_bits_per_dimension`| `int`                | `10`       | Quantization bits per axis for Morton codes               |
| `random_seed`                  | `int`                 | `42`       | RNG seed for reproducibility                             |

## Project Structure

```
fast-patcher/
├── src/
│   ├── clusterAPI.h / .cpp       Main entry point (Patcher::create_clusters)
│   ├── Lloyd.h / .cpp            BFS kernels, seed init, centre update, split
│   ├── ConnectedComponent.h/.cpp Connected-component detection
│   └── bit_array.h               Compact bit-packed boolean array
├── demo/
│   ├── test_cluster.cpp          Demo executable (mesh → clusters → visualise)
│   └── utils/
│       └── polyscope_visualizer.h  Polyscope rendering helpers
├── cmake/recipes/                FetchContent recipes for dependencies
├── CMakeLists.txt                Build configuration
└── README.md                     This file
```

## License

See the repository for license details.
