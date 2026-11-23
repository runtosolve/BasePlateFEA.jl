# BasePlateFEA.jl

Shell finite element model of a cold-formed steel rack upright on a base plate anchored to a
concrete block, for the moment-rotation response of the column base. It is the Julia
implementation of S. Ádány's MATLAB column-base model (`bp_Barcelona.m`, `to_Cris_Nov6_2026`),
extended with metal plasticity and a 3D results viewer. Units are kip and inch.

```julia
using BasePlateFEA

inp = BasePlateInput(geometry = barcelona_geometry(:C),                 # or BasePlateGeometry(W = ..., L = ..., ...)
                     plate_steel = bilinear_steel(53.2, 64.3, 0.15),    # fy, fu, strain at fu; omit for elastic steel
                     upright_steel = bilinear_steel(68.9, 77.8, 0.10))
res = analyze(inp; N = 10.0, Mmax = 100.0, θmax = 0.05)                 # axial preload, then moment
θ, M = moment_rotation(res)
write_viewer("typeC.html", res; title = "Type C, 10 kip")               # 3D viewer page
```

## Model

| Part | Implementation | MATLAB source |
|---|---|---|
| Base plate mesh | DistMesh triangles; upright footprint, hole polygons and hole centres are fixed nodes; elements inside the holes get 1/100 of the thickness | `main_modelgen_dm.m`, `distmesh.m` |
| Upright | structured quadrilaterals on the five walls of the lipped C (centerline dimensions) | `mesh4_rect.m` |
| Shell elements | TriShellFiniteElement.jl and QuadShellFiniteElement.jl. With `Cs = 0` and `drilling = :penalty` (defaults here) these are the MATLAB t31/q42 elements | `stiffmat_e.m` |
| Rigid links | stiff Euler-Bernoulli beams from a master node to the top cross-section, and from each hole centre to the hole perimeter | `ke_beam_EB.m` |
| Concrete | elastic half-space: inverse of the Boussinesq flexibility of the base plate nodes | `stiffmat_concr.m` |
| Contact | compression-only link between every plate node and its concrete dof, elastic to a yield stress then hardening, with elastic unloading | `stifmat_cond.m`, `spr_*.m` |
| Anchor bolts | tension-only spring at each hole centre, same law | `spring_support.m` |
| In-plane restraint | soft springs on the hole perimeter and on every plate node | `fem_support.m` |
| Loading | step 1 axial force; step 2 `V = M/Ht` at height `H` plus the moment `M − V H`, bending about the upright's axis of symmetry | `bp_Barcelona.m` |
| Rotation | `(w₂ − w₁)/h` between two flange points at height `Hm` | `collect_control.m` |

All defaults of `BasePlateInput` are the November 2026 MATLAB values; see its docstring.

### Differences from the MATLAB code

* **Equilibrium iteration.** Every increment is iterated to equilibrium with the consistent
  tangent (residual tolerance 1e-6). The MATLAB code accepts an increment when 0.3 × residual is
  below 0.2 % of the load, usually after one pass.
* **Displacement control.** Step 2 is driven by the work-conjugate displacement of the load
  pattern, so the analysis continues along a plateau or past a peak moment.
* **Plasticity.** Parts with a plastic `Steel` use the geometrically nonlinear, layered von Mises
  elements of the Tri/Quad shell packages (5 Simpson layers by default). The element data are
  built element by element in `src/shell.jl` from those packages' internal functions, because the
  model mixes element types, thicknesses and materials.
* Quirks of the MATLAB code that are reproduced on purpose: the rigid-link torsional stiffness
  `G J / L³`, the contact law scaled by the mean nodal area (`link_area = :average`), the
  Boussinesq patch of 1/10 of √(nodal area), and dropping concrete stiffness entries below 1 % of
  the largest.

## Viewer

`write_viewer(path, result; curves = [...])` writes one HTML page (three.js from a CDN): deformed
model, contours of von Mises stress (top, bottom, mid-surface), plastic strain, membrane and
surface stresses, contact stress on the plate and vertical displacement; anchor bolt forces; and
the moment-rotation and axial force-shortening curves with a load step slider. Test or reference
curves can be overlaid.

## Validation

`validation/` holds the scripts; data and results live outside the repository.

* `compare_octave.jl`: on the meshes exported from the MATLAB code run in GNU Octave (types A, B,
  C at 5, 10, 15, 20 kip), the shell, beam, support and concrete stiffness matrices agree to
  round-off and the initial stiffness and axial shortening to 5 digits. After uplift the Julia
  rotations are larger; one Newton correction of the MATLAB end state reproduces the Julia result.
* `run_barcelona.jl`, `make_report.jl`, `make_viewers.jl`: Barcelona tests (Roure et al.), types
  A, B, C, with measured thicknesses and yield/tensile strengths.

## Layout

```
src/input.jl     geometry, materials, options      src/shell.jl    shell elements (linear and plastic)
src/distmesh.jl  DistMesh                           src/solve.jl    model assembly and analysis
src/mesh.jl      model generation                   src/results.jl  element stresses
src/parts.jl     links, supports, concrete, springs src/viewer.jl   viewer export (viewer/template.html)
archive/         the June 2026 single-file port and its scripts
```
