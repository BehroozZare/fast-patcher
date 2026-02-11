## Fast Graph Patching (Plan)

### Goal

Partition a graph into **patches (clusters)** such that each patch has size **≤ `patch_size`**, using a Lloyd-style iterative refinement (k-means analog on graphs).

---

## Inputs

### Required

* **Graph** `G` in **CSR** format

  * Undirected graph, but **CSR is provided without diagonal/self edges**.
* **`patch_size`**: maximum desired patch size.

### Optional

* **`opt`**: configuration struct (can be `nullptr` for now). Planned fields:

  * `int lloyd_iters_to_add_seed = 5;`
  * `int random_seed = ...;` (for reproducible “random pick” in center update)
  * (later: BFS backend, queue type, parallelism options, etc.)

---

## Outputs

* `vertex_to_patch[v]`: patch id for each vertex `v` (patch ids can be globally unique across CCs).

---

## Algorithm Steps

### Step 0 — Load inputs

**Input:** CSR graph (no diagonals), `patch_size`, optional `opt`
**Output:** internal graph representation + config defaults

---

### Step 1 — Connected components (CC)

Compute connected components to avoid clustering across disconnected subgraphs.

**Output:**

* `vertex_to_cc[v]`: CC id for each vertex

---

### Step 2 — Choose number of seeds per CC

For each connected component `cc`:

* Let `n_cc = |V_cc|`.
* Compute initial number of seeds (clusters), e.g.:

  * `k_cc = ceil(n_cc / patch_size)`
    (This is the “minimum” number of patches needed if all were perfectly balanced.)

**Output:**

* `k_cc` for each CC

---

### Step 3 — Initialize seeds (cluster centers)

For each CC, pick `k_cc` initial seeds (centers). You can start simple:

* random vertices in the CC, or

**Output:**

* `seeds_cc[i]` for `i = 0..k_cc-1`

---

## Step 4 — Lloyd-style patch refinement (per CC)

Run iterative refinement for each connected component independently.

### Step 4.1 — Distance computation (BFS from seeds)

For each seed, compute graph distances to all vertices in the CC using BFS.

Two equivalent implementations:

* **Multi-source BFS with labels:** initialize queue with all seeds and propagate the nearest seed label

**Output:**

* `dist[v]`: distance to assigned seed (or temporary)
* `label[v]`: which seed/cluster currently “owns” `v`

---

### Step 4.2 — Cluster assignment

Assign each vertex `v` to the cluster of the seed with **minimum distance**:

* `label[v] = argmin_i dist_i[v]`

Tie-breaking policy (make deterministic if you want reproducibility):

* for tie, pick the cluster (patch) with smaller size

**Output:**

* `clusters[i]` (implicit via `label[v]`)

---

### Step 4.3 — Update cluster centers (graph “k-means center”)

For each cluster `C`:

1. **Find boundary vertices** of `C`:
   vertices in `C` that have at least one outgoing edge to a vertex in a different cluster.
2. Run a **multi-source BFS starting from the boundary vertices**, restricted to vertices in `C`, in “waves”.
3. Let `W_last` be the last wave reached (vertices farthest from the boundary within the cluster).
4. Pick the new center as:

   * a **random vertex from `W_last`** (random but controlled by `opt.random_seed`), or
   * optionally a deterministic choice (min id) if you prefer.

**Output:**

* updated `seeds_cc[i]` for next iteration

---

### Step 4.4 — Handle oversized patches (optional split)

After a fixed number of Lloyd iterations (`opt.lloyd_iters`, default 5), check patch sizes:

* If any cluster has `|C| > patch_size`, then **add an extra seed inside that cluster** and continue refinement.

  * New seed selection can reuse the “deep interior” concept:

    * pick from `W_last` (farthest-from-boundary wave) of that oversized cluster
* Continue until the convergence condition holds:

**Convergence condition:**
All clusters satisfy `|C| ≤ patch_size`.

**Output:**

* final seeds and assignments satisfying the size constraint (when possible)

---

## Step 5 — Final labeling

Combine results across CCs into global patch ids:

**Output:**

* `vertex_to_patch[v]` for all vertices