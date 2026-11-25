# File formats

All quantities are in SI units (m, kg, s, N, Pa, W) unless stated otherwise. Text files are
read line by line; empty lines and lines starting with `#` are skipped, and numbers may be
separated by spaces or commas.

## Parameter file

One `key = value` per line; `#` starts a comment. Unknown keys stop the program with an error.
The same keys can be given on the command line (`tillbridge run ... --set key=value ...`), in a
calibration file (`set key value`) or through `tb::Params` in the library. Every run writes the
parameters it used to `<out>_params.txt`.

| Key | Default | Meaning |
|---|---|---|
| `E` | 2.6e6 | Young's modulus of the particles (Pa) |
| `nu` | 0.3 | Poisson's ratio of the particles |
| `mu_pp` | 0.6 | particle–particle friction coefficient |
| `mu_pw` | 0.5 | particle–wall and particle–tool friction coefficient |
| `e` | 0.4 | coefficient of restitution (sets the damping ratio) |
| `mu_r` | 0.3 | rolling resistance coefficient (particle–particle, type A) |
| `k_coh` | 0 | cohesion energy density (Pa): attractive force k_coh·π·R*·δ |
| `E_stone` | -1 | Young's modulus of material class 1 (stones); ≤ 0: same as `E` |
| `E_wall`, `nu_wall` | 2.1e11, 0.3 | elastic constants of the walls and the tool |
| `g` | 9.81 | gravity (along −z) |
| `dt` | 1e-4 | time step (s) |
| `t_end` | 1.0 | end time (s) |
| `t_tool` | 0.1 | time when the tool starts to move (s) |
| `v_tool` | 1.5 | tool speed along +x (m/s) |
| `out_dt` | 0.01 | output interval of the force history (s) |
| `Lx`, `Ly`, `Lz` | 1, 1, 0.3 | box size; walls at x = 0, Lx, y = 0, Ly and floor z = 0 |
| `avg_x0`, `avg_x1` | – | averaging window: tool displacement between these values (m) |
| `v_max` | 50 | the run stops (exit code 3, `unstable = 1` in the summary) if a particle exceeds this speed (m/s) |
| `ke_growth` | 0 | if > 0: the run also stops as unstable when, before the tool moves, the kinetic energy exceeds `ke_growth` times its minimum since t = 0.02 s plus 0.005 J/kg of bed mass (a resting bed must not gain energy; for cohesive beds, not for settling runs) |
| `snap_dt` | 0 | if > 0: particle positions and velocities every `snap_dt` seconds to `<out>_snap_NNNN.txt` |
| `bond_E` | 0 | Young's modulus of the parallel bonds (Pa); 0 switches bonds off |
| `bond_kr` | 2.5 | normal to shear stiffness ratio of the bonds |
| `bond_lam` | 0.5 | bond radius as a fraction of the smaller sphere radius |
| `bond_sc`, `bond_tc` | 1e30 | tensile and shear strength of the bonds (Pa) |
| `bond_gap` | 0.05 | bonds are created at t = 0 between spheres closer than (r_i + r_j)(1 + `bond_gap`) |
| `bond_zeta` | 0.2 | damping ratio of the bonds |
| `bond_zmin` | −1e9 | no bonds to bodies whose centre is below this height (m) |
| `out` | tbdem | output prefix (set by the command line program from its OUT argument) |

The time step should be below about 20 % of the Rayleigh time
T_R = π R (ρ/G)^0.5 / (0.1631 ν + 0.8766) of the smallest particle.

## Bed file

```
NB NS
x y z m I ns [mat [vx vy vz [wx wy wz]]]
dx dy dz r            (ns lines)
...
```

Body centre, mass, moment of inertia (spherical), number of spheres, material class
(0 soil, 1 stone, 2 fixed), optional initial velocity and angular velocity; then the sphere
centres relative to the body centre and their radii. NS is the total number of spheres. Bodies
of class 2 do not move (anchored subsoil, see `tillbridge prepare --fix-below`).

## Tool file

```
NT
ax ay az bx by bz cx cy cz      (NT lines)
```

Triangles in the tool frame. The tool frame moves with the tool: at time t ≥ t_tool its origin
is at (v_tool·(t − t_tool), 0, 0) in the box frame. Contact is two-sided. `tillbridge tool`
converts ASCII or binary STL files (with a unit factor, a translation and optional refinement).
An empty tool (`0`, or the word `none` on the command line) is used for settling runs.

## Outputs of a run

| File | Content |
|---|---|
| `<out>_force.csv` | t, Fx, Fy, Fz (force of the soil on the tool), kinetic energy of the particles (translational and rotational, J), tool displacement |
| `<out>_tri.csv` | per triangle, means over the averaging window: normal contact force (positive in compression), number of contacts, force vector of the soil on the triangle |
| `<out>_triW.csv` | per triangle, means over the averaging window: power dissipated by the tangential (friction) and by the normal (damping) contact forces, W |
| `<out>_state.txt` | final bed in the bed format (sphere offsets rotated) |
| `<out>_vel.txt` | final velocity and angular velocity of each body |
| `<out>_depart.txt` | particles leaving the tool: sphere index, body index, position, velocity and tool displacement at the end of the last tool contact |
| `<out>_bonds.txt` | with bonds: body indices of the bonded pair and 1 (intact) or 0 (broken) |
| `<out>_snap_NNNN.txt` | with `snap_dt` > 0: body positions and velocities |
| `<out>_summary.txt` | `key = value`: version, run size, steps, CPU time, mean number of contacts, bonds created and broken, number of samples in the averaging window and the mean tool force; since version 2.1 also the energy accounting over the averaging window (`E_distance` in m; `E_tool`, `E_tool_contacts`, `E_bulk`, `dE_kinetic`, `dE_potential`, `E_closure` in J, see docs/model.md); for an unstable run the time and the speed or energy at which it stopped |
| `<out>_params.txt` | the parameters of the run |

The files have the same layout as those of the earlier C solver (version 1.x), so results of
both can be post-processed with the same commands.

## Post-processing outputs

* `tillbridge wear PREFIX TOOL --out FILE.csv`: per triangle the mean normal pressure (Pa),
  frictional power density q_f (W/m²) and mean slip speed (m/s, for triangles with a pressure
  above `--p-min`); the totals are printed. q_f is a relative wear index: by Archard's law the
  wear rate is K p v_s / H = K q_f / (μ H).
* `tillbridge energy PREFIX TOOL --speed V [--bed INITIAL --x0 X0 --x1 X1]`: work per metre of
  travel (J/m) — draft, friction on the tool, normal damping, the rest (soil and slice) and,
  with the initial bed, kinetic energy and lifting work of the soil leaving the tool.
* `tillbridge fem ... --vtk FILE`: displacements and von Mises stress at the nodes (stress
  components extrapolated from the Gauss points to the corner nodes, mid-side nodes as the mean
  of their edge, averaged over the elements at a node, as in CalculiX) and the largest von Mises
  stress at the Gauss points of every element, in legacy VTK format for ParaView.
* `tillbridge fem ... --ccx FILE.inp`: the same model for CalculiX (`*NODE`, `*ELEMENT` with
  C3D8, C3D20R or C3D10, node sets `FIXED` and `LOADED`, `*ELASTIC`, the nodal forces of the
  distributed coupling as `*CLOAD`; nodes numbered from 1 in the order of the solver). Run it
  with `ccx -i FILE` and compare the `.frd` results with the VTK file.

## Gmsh meshes

`tillbridge fem`, `tillbridge convert` and `tb::read_gmsh` read ASCII `.msh` files of versions
2.2 and 4.1 with 8-node and 20-node hexahedra and 10-node tetrahedra (Gmsh element types 5,
17 and 11), converted to the node order of the solver. Every physical group, of any dimension,
becomes a node component named after the group (unnamed groups: `PHYS<dim>_<tag>`), so the
clamped face and the loaded face are marked as physical surfaces in Gmsh. The file holds no
material, so `--E` (and `--nu`, default 0.3) are given on the command line.
`tillbridge convert MODEL OUT.msh` writes version 2.2 with the solid elements in the physical
volume `SOLID` and every node component as a physical group of one-node point elements.

## ANSYS meshes

`tillbridge fem` and `tb::AnsysModel` read node blocks (`NBLOCK`), solid element blocks
(`EBLOCK`) with 8-, 20- and 10-node elements (SOLID185, SOLID186 and SOLID187 node layouts;
20-node elements are integrated with the reduced 2×2×2 rule), node components (`CMBLOCK`), `MP`
material lines and scalar parameters, from archive files (`.cdb`, `CDBWRITE` blocked format)
and Mechanical input files (`.dat`). Units are those of the model; the
command line program applies forces in N and moments in N·mm for a model in mm.
