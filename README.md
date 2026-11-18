# TillBridge

Calibrated discrete element loads on tillage tools, carried to finite element stresses and fatigue life.

TillBridge links three steps that are usually done in separate programs when a tillage tool is
designed from simulation:

1. **Discrete element simulation of the soil–tool interaction.** Rigid multi-sphere particles
   (spheres or volume-equivalent clods) in a soil bin are loaded by a rigid triangulated tool
   moving at constant speed. Contacts follow the Hertz–Mindlin model with damping from the
   coefficient of restitution, Coulomb friction with tangential history, linear cohesion,
   type-A rolling resistance and, optionally, parallel bonds. Velocities are evaluated at the
   contact point, so particle rotation and sliding are coupled. Per-triangle normal and
   tangential forces and contact powers are recorded on the tool.
2. **Calibration of the soil parameters.** Gaussian-process surrogates and Bayesian
   optimisation fit the DEM parameters to measured quantities (target matching) or minimise
   one quantity while others stay in given intervals (constrained calibration). Runs are
   proposed in batches and logged, so a calibration can be resumed.
3. **Assessment of the tool.** The work of the draft force is split into friction on the tool,
   normal contact damping and dissipation in the soil; the frictional power density on each
   triangle gives a wear index (Archard). The soil load is reduced to a wrench at the centre of
   pressure, transferred to the finite element model of the shank by a distributed coupling,
   and solved with 8/20-node hexahedra or 10-node tetrahedra and an algebraic-multigrid
   preconditioned conjugate gradient solver. The force history gives a stress history, counted
   with the rainflow method and summed with the Palmgren–Miner rule on the EN 1993-1-9 S–N curves.

The code is C++17 without external dependencies; OpenMP is used when the compiler provides it.
No commercial DEM or FE program is needed. Meshes are read from the open Gmsh format (`.msh`) or
from ANSYS (`.cdb`, Mechanical `.dat`), STL tool geometry from any CAD program. `tillbridge convert`
passes an ANSYS model on to Gmsh, and `tillbridge fem --ccx` writes the same model as a CalculiX
input file for an independent check of the solution.

## Build

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

Requirements: CMake ≥ 3.14 and a C++17 compiler (GCC ≥ 8, Clang ≥ 7, MSVC 2019 or newer).
This builds the command line program `tillbridge`, the static library `libtillbridge` with its
headers in `include/tillbridge/`, the worked example `chisel_example` and the tests.
Options: `-DTILLBRIDGE_OPENMP=OFF`, `-DTILLBRIDGE_TESTS=OFF`, `-DTILLBRIDGE_EXAMPLES=OFF`.

## Quick start

```
tillbridge bed column.txt --Lx 1.2 --Ly 0.6 --depth 0.25 --radius 0.015
tillbridge run column.txt none settle --set t_end=1.5 t_tool=99 Lx=1.2 Ly=0.6 Lz=0.3
tillbridge prepare settle_state.txt bed.txt --depth 0.25
tillbridge tool tine.stl tine.txt --unit 0.001 --translate -0.1 0.3 0
tillbridge run bed.txt tine.txt run1 --set v_tool=1.5 t_tool=0.1 t_end=0.8 Lx=1.2 Ly=0.6 avg_x0=0.3 avg_x1=0.9
tillbridge energy run1 tine.txt --speed 1.5
tillbridge wear run1 tine.txt --smooth 1 --out run1_wear.csv
tillbridge loads run1 tine.txt --ref 0 0.3 0.45
tillbridge convert shank.cdb shank.msh
tillbridge fem shank.msh --E 200000 --fixed CLAMP --loaded FACE --force -850 40 310 --moment 0 1.2e5 0 --vtk shank.vtk
tillbridge fatigue run1 --x0 0.3 --x1 0.9 --sigma-ref 85 --detail 71 --speed 1.5
tillbridge calibrate calibration.txt
```

`tillbridge help` lists all commands and options. Results are printed as JSON-like key–value
lists; the solver writes text and CSV files described in `docs/file_formats.md`. The
calibration file is described in `docs/calibration.md`.

The library can be used directly; `examples/chisel/chisel_example.cpp` goes from a generated
soil bed through a synthetic calibration with a known answer to the FE model of the tine and its
fatigue life. Python is used only for the figures (`examples/chisel/make_figures.py`, NumPy
and Matplotlib).

## Verification

`ctest` runs the unit tests (`tillbridge_tests`, about 15 s) and a run of every command of the
command line program:

* contact-law benchmarks with analytical solutions (after Chung and Ooi 2011): static Hertz
  overlap, normal restitution and its independence of the impact speed, oblique impacts on a
  wall and between particles in the sliding regime, a spinning sphere set down on a floor
  (v = 2/7 R ω₀), a sphere dragged by a moving tool plate (v = 2/7 v_tool), absence of rolling
  resistance on walls, parallel bonds that hold or break, closest-point geometry;
* regression against the reference C solver from which the code was ported (chisel tine in a
  settled bed of 622 spheres, `tests/data/regression`): with GCC or Clang and one thread the
  force history is reproduced to all printed digits; otherwise the mean forces must agree
  within 5 %;
* finite elements: patch tests (distorted hex20, tet10), cantilever deflection against
  Timoshenko beam theory, bending stress, AMG against Jacobi preconditioning, ANSYS and Gmsh readers;
* fatigue: rainflow counts of the ASTM E1049 example, S–N curve knee points (2·10⁶, 5·10⁶ and
  10⁸ cycles), Miner scaling;
* load transfer (wrench and distributed coupling), bed and tool input/output (text, ASCII and
  binary STL), clod volume, parameter files, Gaussian process and expected improvement, and
  recovery of known parameters by both calibration modes and after a restart from the log.

Continuous integration builds and tests on Linux (GCC, Clang), macOS (Apple Clang) and Windows (MSVC).

## Files

| Path | Content |
|---|---|
| `include/tillbridge/` | library headers: `dem`, `bed`, `tool`, `params`, `calibrate`, `post` (energy, wear, load transfer, fatigue), `fem`, `sparse`, `ansys`, `gmsh`, `vec` |
| `src/` | library sources |
| `app/main.cpp` | command line program |
| `docs/` | file formats, calibration file, model description |
| `examples/chisel/` | worked example and figure script |
| `examples/plough/` | mouldboard plough of the paper: soil bin, tool, shank mesh (Gmsh), commands and the check against CalculiX |
| `examples/benchmarks/` | verification cases, FE convergence and solver scaling |
| `tests/` | verification tests and test data |

## Citation

See `CITATION.cff`.

## Licence

MIT — see `LICENSE.txt`.
