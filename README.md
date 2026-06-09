# Column Base Plate FEA Package

Julia package for nonlinear incremental finite element analysis of thin-walled cold-formed steel column base plate connections.  
Converted from MATLAB (S. Adany, 2026).

---

## Table of Contents

- [Overview](#overview)
- [Directory Structure](#directory-structure)
- [How to Run](#how-to-run)
- [Dependencies](#dependencies)
- [Key Data Structures](#key-data-structures)
- [Coordinate Systems](#coordinate-systems)
- [Source Files](#source-files)
  - [Types](#types)
  - [Cross-Section Properties](#cross-section-properties)
  - [Finite Element Properties](#finite-element-properties)
  - [Coordinate Transformations](#coordinate-transformations)
  - [Gauss Integration](#gauss-integration)
  - [Element Stiffness Matrices — In-Plane (Membrane)](#element-stiffness-matrices--in-plane-membrane)
  - [Element Stiffness Matrices — Out-of-Plane (Bending)](#element-stiffness-matrices--out-of-plane-bending)
  - [Stress Calculation](#stress-calculation)
  - [Load Vectors](#load-vectors)
  - [Mesh Generation](#mesh-generation)
  - [Global Assembly](#global-assembly)
  - [Support Conditions](#support-conditions)
  - [Main Analysis Functions](#main-analysis-functions)
  - [Visualization](#visualization)
- [Output Files](#output-files)
- [Analysis Workflow](#analysis-workflow)

---

## Overview

Models a **C-section cold-formed steel column** bolted to a **steel base plate**, which rests on a rigid foundation. The model captures:

- Thin-walled shell behaviour (membrane + bending) using Mindlin plate theory
- Bolt hole geometry approximated as polygons
- One-way (compression-only) conditional springs under the base plate — modelling contact with the foundation
- Bolt bearing modelled as radial springs around the hole perimeter
- Incremental load application to track nonlinear contact response
- Minor-axis bending, major-axis bending, axial force, and shear loads

---

## Directory Structure

```
julia files/
├── Project.toml          # Package definition and dependencies
├── Manifest.toml         # Pinned dependency versions (auto-generated)
├── README.md             # This file
│
├── src/
│   ├── BasePlateFEA.jl        # Main module — imports all source files
│   ├── types.jl          # Struct definitions shared across the package
│   │
│   │── Cross-section
│   ├── cs_def.jl
│   ├── cs_prop.jl
│   ├── cs_data.jl
│   │
│   │── FE element data
│   ├── FE_props.jl
│   │
│   │── Coordinate transforms
│   ├── plate_transform.jl
│   ├── ct_3node_g2e.jl
│   ├── ct_4node_g2e.jl
│   ├── rotate3d.jl
│   ├── add_drill.jl
│   │
│   │── Gauss integration
│   ├── GL3.jl
│   ├── GL4.jl
│   │
│   │── Element stiffness — membrane
│   ├── ke_uv_3n_nocondens_6dof.jl
│   ├── ke_uv_3n_condens_from_8to6dof.jl
│   ├── ke_uv_4n_nocondens_8dof_num.jl
│   ├── ke_uv_4n_condens_from_12to8dof_num.jl
│   │
│   │── Element stiffness — bending (Mindlin)
│   ├── ke_wt_3n_condens_from_12to9dof_MIN.jl
│   ├── ke_wt_4n_condens_from_18to12dof_num.jl
│   │
│   │── Stress output
│   ├── stress_uv_3n.jl
│   ├── stress_uv_4n.jl
│   │
│   │── Load vectors
│   ├── lv_el3_edge.jl
│   ├── lv_el4_edge.jl
│   ├── load_vec_UDL.jl
│   ├── load_vec_edge.jl
│   ├── load_vec_point.jl
│   │
│   │── Mesh generation
│   ├── collect_anchor.jl
│   ├── mesh4_rect.jl
│   ├── edges_by_plate.jl
│   ├── elem_divide_by_plate.jl
│   ├── four2three_by_plate.jl
│   │
│   │── Global assembly
│   ├── stiffmat_e.jl
│   ├── fem_support.jl
│   ├── cond_support.jl
│   │
│   │── Top-level analysis
│   ├── main_modelgen.jl
│   ├── main_stressout.jl
│   ├── FEanal_incr.jl
│   │
│   └── Visualization
│       └── plot_res_tri.jl
│
└── scripts/
    └── bp_test.jl        # Main analysis script (C-section base plate example)
```

---

## How to Run

```powershell
# from the julia files/ directory
julia --project=. scripts/bp_test.jl
```

Or interactively in the REPL:

```julia
using Pkg
Pkg.activate(".")          # activate the project environment
include("scripts/bp_test.jl")
```

---

## Dependencies

| Package | Role |
|---|---|
| `LinearAlgebra` | Matrix operations, `cross`, `norm`, `det`, `inv`, `diagm` |
| `SparseArrays` | Sparse stiffness matrix assembly |
| `Statistics` | `mean` for node averaging |
| `Random` | Random quad→triangle splitting |
| `GLMakie` | 3D mesh visualization with filled color patches |
| `GeometryBasics` | `Point3f`, `TriangleFace`, `Mesh` types for GLMakie |

All dependencies are pinned in `Manifest.toml` and installed automatically via `Pkg.instantiate()`.

---

## Key Data Structures

### `Tolerances`
```julia
mutable struct Tolerances
    zero_thickn  :: Float64   # minimum thickness (in)
    zero_spring  :: Float64   # minimum spring stiffness (kip/in)
    zero_dist    :: Float64   # node coincidence tolerance (in)
    zeroFElength :: Float64   # minimum element size (in)
end
```

### `FEData`
Describes one finite element type (stored in a `Dict{Int, FEData}` keyed by element ID 11, 12, 21, 22).
```julia
mutable struct FEData
    shape        :: Int     # 3 = triangle, 4 = quad
    nodenr       :: Int     # nodes per element
    m_dof_ini    :: Int     # membrane DOFs before condensation
    m_dof_final  :: Int     # membrane DOFs after condensation
    m_ifcondens  :: Bool    # membrane: use static condensation?
    m_numint     :: Bool    # membrane: numerical integration?
    m_nGpe       :: Int     # membrane: number of Gauss points
    b_theory     :: String  # bending theory ("Mindlin")
    b_dof_ini    :: Int     # bending DOFs before condensation
    b_dof_final  :: Int     # bending DOFs after condensation
    b_ifcondens  :: Bool
    b_numint     :: Bool
    b_nGpe       :: Int
    nonzer_ke    :: Int     # estimated non-zeros in element Ke
    nonzer_lv    :: Int     # estimated non-zeros in element load vector
    nGpg         :: Int     # Gauss points for global assembly
end
```

### Element type IDs
| ID | Shape | Membrane | Bending |
|---|---|---|---|
| `11` | Triangle 3-node | 6 DOF, analytic, no condensation | Mindlin, 12→9 DOF |
| `12` | Triangle 3-node | 8→6 DOF, analytic, condensed | Mindlin, 12→9 DOF |
| `21` | Quad 4-node | 8 DOF, numeric, no condensation | Mindlin, 18→12 DOF |
| `22` | Quad 4-node | 12→8 DOF, numeric, condensed | Mindlin, 18→12 DOF |

### `HoleData`
Stores polygon approximation of circular bolt holes for one plate.
```julia
mutable struct HoleData
    coord  :: Vector{Vector{Float64}}   # [x_centre, y_centre] per hole
    radius :: Vector{Float64}           # hole radius per hole
    sector :: Vector{Int}               # polygon side count per hole
    holes  :: Vector{Matrix{Float64}}   # polygon vertices (n+1 × 2) per hole
end
```

### Input tables (all `Vector{Vector{Any}}` or `Vector{Vector{Float64}}`)

| Table | Format | Description |
|---|---|---|
| `PLATE_table` | `[x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4 t imat]` | Corner coordinates, thickness, material per plate |
| `MAT_table` | `[Ex Ey nuxy nuyx G rho]` | Orthotropic material constants |
| `LOAD_table` | mixed `Any` entries | `"edge"`, `"UDL"`, or `"point"` load definitions |
| `SUPPORT_table` | mixed `Any` entries | Ordinary spring supports at points, lines, or surfaces |
| `CONDSUP_table` | mixed `Any` entries | Conditional (one-way) spring supports |
| `HOLE_table` | `Vector{Union{Nothing, HoleData}}` | Bolt hole data per plate |

### Node and element arrays

| Array | Size | Columns |
|---|---|---|
| `nodes` | `total_nodes × 5` | `x, y, z, id_orig, id_active` |
| `el_nodes` | `total_elems × 9` | global node indices (up to 9; unused = 0) |
| `el_props` | `total_elems × 5` | `plate_id, fe_type_id, thickness, mat_id, hole_flag` |

- `id_orig > 0` → node is active (canonical); `id_orig = 0` → merged duplicate  
- `id_active` → sequential active node number 1…`nodenr` (used for DOF indexing)  
- `hole_flag = 1` → normal element; `-1` → element inside hole (thickness ÷100, excluded from display)

### DOF ordering
Each active node has `dpn = 6` DOFs ordered: **u, v, w, θ_x, θ_y, θ_z**  
(global X, Y, Z translations + rotations). DOF index for active node `n`, component `k`:  
```
dof = (n - 1) * dpn + k
```

### `out_res` array (post-processing)
Shape: `(nelems, 9, 10, 3)` — element × local node × component × layer (bottom/mid/top)

| Component index | Quantity |
|---|---|
| 1–3 | Global displacements U, V, W |
| 4–6 | Surface-local displacements u, v, w |
| 7 | In-plane stress σ_x |
| 8 | In-plane stress σ_y |
| 9 | In-plane shear stress τ_xy |
| 10 | Conditional support status (contact) |

---

## Coordinate Systems

| Label | Description |
|---|---|
| `g` | Global Cartesian (X, Y, Z) — world frame |
| `s` | Surface-fixed Cartesian — aligned to each plate's plane |
| `e` | Element-local Cartesian — aligned to each finite element |
| `l` | Local parametric (ξ, η) — reference element coordinates |

Transformation matrices:
- `plate_transform` → rotation matrix **g ↔ s**
- `ct_3node_g2e` / `ct_4node_g2e` → rotation matrix **g ↔ e**
- `rotate3d` → block-diagonal rotation for full element DOF vector

---

## Source Files

### Types

#### `types.jl`
Defines the three shared structs used throughout the package:
- `Tolerances` — numerical tolerances (distances, spring thresholds)
- `FEData` — element type descriptor (DOF counts, integration settings)
- `HoleData` — bolt hole polygon geometry

---

### Cross-Section Properties

#### `cs_def.jl` — `cs_def(h, b, c, tb, th, tc)`
Defines node and element connectivity for a thin-walled C-section.  
Returns `cs_node` (node coordinates) and `cs_elem` (connectivity + thickness).

#### `cs_prop.jl` — `cs_prop(cs_node, cs_elem)`
Calculates basic cross-section properties from the node/element arrays:
- Area `A`, centroid (`yCG`, `zCG`)
- Moments of inertia `Iy`, `Iz` (centroidal) and `Iyr`, `Izr` (about reference axis)

#### `cs_data.jl` — `cs_data(h, b, c, tb, th, tc)`
Convenience wrapper: calls `cs_def` then `cs_prop`.  
Returns `A, yCG, zCG, Iyr, Iy, Izr, Iz`.

---

### Finite Element Properties

#### `FE_props.jl` — `FE_props()`
Creates and returns a `Dict{Int, FEData}` with entries for the four supported element types:
- `11` — triangular, no membrane condensation, analytic
- `12` — triangular, membrane condensation 8→6 DOF, analytic
- `21` — quad, no membrane condensation, numeric Gauss
- `22` — quad, membrane condensation 12→8 DOF, numeric Gauss

---

### Coordinate Transformations

#### `plate_transform.jl` — `plate_transform(P1, P2, P3, P4)`
Computes the rotation matrix **T** (3×3) from global to plate-local coordinates for a flat rectangular plate defined by four corner points.  
`T` columns are the three local axes; midpoint `P0` is also returned.

#### `ct_3node_g2e.jl` — `ct_3node_g2e(elems_row, nodes)`
Rotation matrix **T** from global to 3-node triangular element coordinates.  
Also returns element-frame corner coordinates `P1e, P2e, P3e`.

#### `ct_4node_g2e.jl` — `ct_4node_g2e(elems_row, nodes)`
Same as above for a 4-node quadrilateral element.  
Returns `T, P1e, P2e, P3e, P4e`.

#### `rotate3d.jl` — `rotate3d(T3, nodenr, dpn)`
Expands the 3×3 element rotation matrix `T3` into a full block-diagonal transformation matrix `T` of size `(nodenr×dpn) × (nodenr×dpn)`. Used to rotate element stiffness matrices to the global frame.

#### `add_drill.jl` — `add_drill(k5dpn, nodenr, ifelastic)`
Expands a 5-DOF-per-node stiffness matrix to 6-DOF by inserting a drilling DOF (θ_z) at each node. The drilling stiffness is set to 1/100 of the minimum bending stiffness for numerical stability.

---

### Gauss Integration

#### `GL3.jl` — `GL3(np)`
Gauss–Legendre integration data for **triangular** elements.  
Supports `np = 1` (centroid rule) or `np = 4` points.  
Returns integration point coordinates `(locx, locy)` and weights `wei`.

#### `GL4.jl` — `GL4(np)`
Gauss–Legendre integration data for **quadrilateral** elements, domain [−1, +1]².  
Supports `np = 1, 2, 3, 4` (i.e., 1×1 to 4×4 point rules).  
Returns `locx, locy, wei` as `np×np` matrices.

---

### Element Stiffness Matrices — In-Plane (Membrane)

All functions return `ke_uv` — the in-plane stiffness matrix with DOF order `u1 v1 u2 v2 … `.

#### `ke_uv_3n_nocondens_6dof.jl`
3-node triangle, 6 DOF total (2 per node), **analytically integrated**, no condensation.  
Closed-form constant-strain triangle (CST).

#### `ke_uv_3n_condens_from_8to6dof.jl`
3-node triangle, analytically integrated. Starts with 8 DOF (adds a bubble node at centroid), then **statically condenses** to 6 DOF. Improves bending accuracy over plain CST.

#### `ke_uv_4n_nocondens_8dof_num.jl`
4-node quadrilateral, 8 DOF, **numerically integrated** (Gauss), no condensation. Standard isoparametric Q4 element.

#### `ke_uv_4n_condens_from_12to8dof_num.jl`
4-node quad, numerically integrated. Adds 4 mid-side bubble DOFs (12 total), then **statically condenses** to 8 DOF. Provides enhanced membrane performance (similar to Q4 with incompatible modes).

---

### Element Stiffness Matrices — Out-of-Plane (Bending)

All functions return `ke_wt` — the bending stiffness matrix with DOF order `w1 θ_y1 θ_x1 w2 … ` per node.

#### `ke_wt_3n_condens_from_12to9dof_MIN.jl`
3-node triangle, **Mindlin plate theory**, analytically integrated.  
Starts with 12 DOF (adds 3 midside nodes for w), condenses to 9 DOF.  
Includes shear correction. Both shear (`kwts`) and bending (`kwtb`) parts assembled separately.

#### `ke_wt_4n_condens_from_18to12dof_num.jl`
4-node quad, **Mindlin plate theory**, numerically integrated (Gauss).  
Starts with 18 DOF (4 corners × 3 + 2 bubble shape functions), condenses to 12 DOF.  
Shear locking reduced by using extra quadratic shape functions.

---

### Stress Calculation

#### `stress_uv_3n.jl` — `stress_uv_3n(d_elem, P1, P2, P3, t, ...)`
Computes in-plane stresses `[σ_x, σ_y, τ_xy]` for a 3-node element.  
Stress is constant over the element (CST). Returns a `3 × npoint` matrix.

#### `stress_uv_4n.jl` — `stress_uv_4n(d_elem, P1, P2, P3, P4, t, ...)`
Computes in-plane stresses for a 4-node element at corner or Gauss points.  
Uses the isoparametric Jacobian to form the B-matrix at each evaluation point.

---

### Load Vectors

#### `lv_el3_edge.jl` — `lv_el3_edge(elems_row, nodes, px, py, pz, dpn)`
Element load vector for a **3-node triangle** due to linearly varying edge traction.  
`px, py, pz` are load intensities at each of the 6 edge endpoints (2 per edge × 3 edges).

#### `lv_el4_edge.jl` — `lv_el4_edge(elems_row, nodes, px, py, pz, dpn)`
Same for a **4-node quad** (8 edge endpoint values: 2 per edge × 4 edges).

#### `load_vec_UDL.jl` — `load_vec_UDL(...)`
Assembles the global load vector for **uniformly distributed loads (UDL)** over entire plates.  
Integrates analytically using element area integrals.

#### `load_vec_edge.jl` — `load_vec_edge(...)`
Assembles the global load vector for **linearly varying edge loads** along plate boundary lines.  
Iterates over all edges in `pl_edge_nodes`, checks which lie on the load line, interpolates load intensity.

#### `load_vec_point.jl` — `load_vec_point(...)`
Assembles the global load vector for **concentrated point loads and moments** at specific node positions.

---

### Mesh Generation

#### `collect_anchor.jl` — `collect_anchor(PLATE_table, SUPPORT_table, CONDSUP_table, tol)`
Collects **anchor points** (seed points, à la Abaqus) for mesh refinement control.  
Anchor points are added wherever:
- A corner of one plate projects onto an edge of another plate (enforces mesh alignment at T-joints)
- A support or conditional support point lies on a plate edge

Returns `p_anchor` — a vector of matrices (one per plate), each row being an anchor point coordinate.

#### `mesh4_rect.jl` — `mesh4_rect(P1, P2, P3, P4, p_anchor_mat, size_x, size_y, tol)`
Generates a structured **rectangular quad mesh** for a single flat plate.  
Uses anchor points to define sub-rectangles, then subdivides each sub-rectangle uniformly at the target element size.  
Returns `nodes` (N×3) and `elems` (M×4) in local plate coordinates mapped back to global.

#### `edges_by_plate.jl` — `edges_by_plate(pl_elems)`
Builds edge connectivity for each plate:
- `pl_el_edges` — for each element, which edge index (1…edge_tot) each side corresponds to
- `pl_edge_nodes` — for each edge, which two nodes it connects (local + global numbering)

Used by `load_vec_edge` to map edge loads onto elements.

#### `elem_divide_by_plate.jl` — `elem_divide_by_plate(pl_elems, pl_nodes, PLATE_table, HOLE_table, tol)`
Splits quad elements into **two triangles** wherever a bolt hole polygon edge passes through them (diagonal cut). Elements fully inside the hole are flagged but not deleted here.

#### `four2three_by_plate.jl` — `four2three_by_plate(pl_elems, rat_q2t)`
Randomly converts a fraction of quad elements to triangles according to `rat_q2t` (ratio per plate).  
`rat_q2t = 0` → all triangles; `rat_q2t = 1` → all quads (no conversion).

---

### Global Assembly

#### `stiffmat_e.jl` — `stiffmat_e(nodes, el_nodes, el_props, nodenr, MAT_table, FE_dat, mat_sym, dpn)`
Builds the **global elastic stiffness matrix** `Ke` (sparse, `nodenr×dpn` square).

For each element:
1. Selects the right membrane (`ke_uv`) and bending (`ke_wt`) function based on `FE_dat`
2. Assembles 5-DOF-per-node element stiffness from membrane + bending parts
3. Optionally adds drilling DOF to make it 6-DOF-per-node
4. Rotates to global frame via `rotate3d`
5. Scatters into the global sparse matrix using active node DOF indices

#### `fem_support.jl` — `fem_support(Ke, PLATE_table, SUPPORT_table, ...)`
Adds **ordinary spring supports** to the stiffness matrix.  
Supports can be:
- `"point"` — spring at a single node
- `"p-line"` — springs along a line
- `"p-surf"` — springs over a rectangular area
- `"surf"` — distributed surface springs (area-weighted, shared between element nodes)

Rigid supports are applied as very large springs (`max(Ke) × 1000`). Returns sparse `Ke_sup`.

#### `cond_support.jl` — `cond_support(Ke, ..., d, dpn, uncond, tolers)`
Adds **conditional (one-way) spring supports** — active only when displacement is in the specified direction (compression-only for contact, tension-only for bolt bearing).  
- On first call (`uncond=true`): applies a small fraction of the spring stiffness unconditionally (to get a starting displacement estimate)
- Subsequently (`uncond=false`): checks each node's displacement `d` to decide whether each spring is active

---

### Main Analysis Functions

#### `main_modelgen.jl` — `main_modelgen(...)`
**Top-level mesh generation.** Orchestrates the full model build:
1. Collects anchor points
2. Meshes each plate with `mesh4_rect`
3. Divides elements at hole boundaries
4. Converts some quads to triangles
5. Builds edge connectivity
6. Merges all plates into global node/element arrays
7. Marks hole elements and reduces their thickness
8. Identifies and deactivates duplicate nodes at plate intersections
9. Assigns sequential active-node IDs

Returns: `nodes, el_nodes, el_props, nodenr, pl_el_edges, pl_edge_nodes`

#### `FEanal_incr.jl` — `FEanal_incr(iincr, dincr, d_ini, Ke, lv, ..., tolers)`
**Incremental FE solver.** On the first increment (`iincr=1`) only:
1. Assembles `Ke` via `stiffmat_e`
2. Adds ordinary supports via `fem_support`
3. Assembles load vector via `load_vec_UDL` + `load_vec_edge` + `load_vec_point`

Every increment:
1. Scales load vector: `dlv = lv × dincr[iincr]`
2. Evaluates conditional supports based on current displacement `d_ini`
3. Solves: `dd = (Ke + Ke_cond) \ dlv`
4. Accumulates: `d_stat = d_ini + dd`

Returns updated `d_stat, Ke, lv`.

#### `main_stressout.jl` — `main_stressout(d, ..., comp, mat_sym, dpn)`
**Post-processing.** Loops over all elements and computes:
- **Displacements** (comp 1–6): global and surface-local at each element node
- **Stresses** (comp 7–9): in-plane σ_x, σ_y, τ_xy, rotated from element frame to plate surface frame
- **Contact status** (comp 10): which nodes have active conditional supports

Returns `out_res` array of shape `(nelems, 9, 10, 3)`.

---

### Visualization

#### `plot_res_tri.jl` — `plot_res_tri(el_nodes, el_props, out_grid, FE_dat, out_res, scale, layer, comp, gridon; ...)`
Renders the FE mesh as **filled, color-mapped 3D patches** using GLMakie.

**Algorithm:**
1. Single pass over all elements → build triangulated face list (quads split 1→3, 1→3→4) + accumulate per-node values
2. Average result values at shared nodes
3. Displace node positions: `def_grid = positions + scale × displacements`
4. Build `GeometryBasics.Mesh` from `Point3f` vertices and `TriangleFace` connectivity
5. Render with `mesh!` (smooth per-vertex color interpolation) + optional `wireframe!` overlay + `Colorbar`

**Key arguments:**

| Argument | Description |
|---|---|
| `scale` | Displacement magnification factor (0 = undeformed) |
| `layer` | 1 = bottom, 2 = middle, 3 = top surface |
| `comp` | Result component index (3=W, 8=σ_y, 10=contact, etc.) |
| `gridon` | `true` to show grey mesh edge lines |
| `title_str` | Plot title and colorbar label |
| `colormap` | Makie colormap symbol (default `:jet`) |
| `azimuth` | Camera azimuth angle in radians |
| `elevation` | Camera elevation angle in radians |

---

## Output Files

All files are saved in the `julia files/` root directory when running `bp_test.jl`:

| File | Description |
|---|---|
| `displacement_w.png` | 3D mesh plot: vertical displacement W (comp=3), deformed shape |
| `stress_sy.png` | 3D mesh plot: in-plane stress σ_y (comp=8), deformed shape |
| `contact.png` | 3D mesh plot: contact status under baseplate (comp=10), undeformed |
| `force_disp.png` | 2D line plot: moment vs. max baseplate displacement |
| `force_disp_data.csv` | Raw force-displacement data (load factor, min/max disp, moment) |

---

## Analysis Workflow

```
bp_test.jl
│
├── Define geometry, materials, loads, holes, supports
│
├── main_modelgen(...)
│   ├── collect_anchor        → seed points for mesh alignment
│   ├── mesh4_rect            → structured quad mesh per plate
│   ├── elem_divide_by_plate  → split quads at hole polygon edges
│   ├── four2three_by_plate   → random quad→triangle conversion
│   ├── edges_by_plate        → edge connectivity table
│   └── deduplicate nodes     → merge coincident plate-junction nodes
│
├── for iincr = 1 … nincr
│   └── FEanal_incr(...)
│       ├── [iincr=1] stiffmat_e      → global elastic K (sparse)
│       ├── [iincr=1] fem_support     → add ordinary spring DOFs to K
│       ├── [iincr=1] load_vec_*      → assemble global load vector
│       ├── cond_support              → activate one-way springs based on d
│       └── solve (K + K_cond) d = f  → displacement increment
│
├── main_stressout(...)
│   ├── Compute displacements at element nodes
│   ├── Compute in-plane stresses (rotate e→s frame)
│   └── Evaluate contact status from conditional support activity
│
└── plot_res_tri(...)  ×3  +  force-displacement chart
    └── Save PNG files
```
