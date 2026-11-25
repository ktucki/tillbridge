# Changelog

## 2.1.0 (2026-11-25)

* Energy accounting in every DEM run: work of the tool, energy taken up at the tool contacts and
  by all other contacts (particle–particle, walls, rolling resistance, bonds), change of kinetic
  and potential energy, and the closure error over the averaging window (`_summary.txt`,
  `EnergyAccount`); `tillbridge energy` prints the directly accounted work done in the soil beside
  the remainder of the draft work, and the closure error (docs/model.md).
* `facet_wrench`: exact wrench of the facet forces in the frame of an FE model (Eq. (8) of the
  paper), used by the chisel example instead of the resultant at the centre of pressure.
* `paper/`: commands that reproduce the examples and figures of the paper, expected results with
  tolerances, and the data of all figures and tables as CSV (`paper/figure_data`).
* Tests: energy closure of the regression case, facet wrench.
* Documentation: one-way coupling and load transfer, fatigue block and life, wear index
  (docs/model.md).

## 2.0.1 (2026-11-18)

* Missing `#include <algorithm>` in the command line program and two test files; the code now
  builds with GCC 13 and later and with MSVC.

## 2.0.0 (2026-11-18)

* Open mesh format: Gmsh `.msh` (versions 2.2 and 4.1) is read by `fem`, and `convert` writes an
  ANSYS model as `.msh`; `fem --ccx` writes the same model as a CalculiX input file, used to
  check the shank of the plough example against an independent solver.
* Nodal stresses are extrapolated from the Gauss points to the corner nodes and averaged, as in
  CalculiX and ANSYS; `--vtk` also writes the largest von Mises stress of every element.
* `examples/plough`: the mouldboard plough of the paper with all inputs and commands.

* The whole program is written in C++17 as one library with a command line program; no
  Python or other dependencies are needed for the computations.
* DEM solver ported from the C solver of version 1.3; with the same compiler and one thread it
  reproduces the 1.3 results to all printed digits (regression test in `tests/data/regression`).
  Output files keep the 1.x layout.
* Parallel bonds (Potyondy and Cundall 2004), anchored subsoil (material class 2), particle
  snapshots, per-triangle contact powers (`_triW.csv`) and departure records of the soil leaving
  the tool (`_depart.txt`).
* New commands `energy` (energy balance of the draft work) and `wear` (pressure, slip speed and
  frictional power density on the tool); `shape` makes parametric chisel and mouldboard bodies.
* Finite elements, algebraic multigrid, ANSYS reader, rainflow and EN 1993-1-9 fatigue, and
  Gaussian-process calibration with a restartable text log, all in C++.
* Calibration for scattered DEM results: one Gaussian process with a noise term per measured
  quantity, expected improvement of the expected misfit, an estimate from the surrogate means
  (`Calibration::estimate`) and surrogate predictions for maps of the consistent region
  (`Calibration::predict`).
* Unit tests with analytical solutions and a command line test run by `ctest`; CI on Linux,
  macOS and Windows.

## 1.0.0 (2026-11-11, not released)

* C solver `tbdem` with a Python package for beds, tools, calibration, load transfer, finite
  elements and fatigue.
