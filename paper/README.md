# Reproducing the paper

Commands, expected results and data of the paper
*TillBridge: an open discrete element–finite element chain for tillage tools, from calibrated
soil loads to energy, wear and fatigue life* (SoftwareX), for release v2.1.0.

## Build and test

```
git clone --branch v2.1.0 https://github.com/ktucki/tillbridge
cmake -S tillbridge -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

`ctest` runs the unit tests (`tillbridge_tests`) and a run of every command of the program.
All tests must pass on Linux, macOS and Windows.

## Chisel example and Fig. 7 (Section 4.2)

```
export OMP_NUM_THREADS=2            # thread count used in the paper (Windows: set OMP_NUM_THREADS=2)
build/chisel_example out            # bed, four reference passes, calibration, FE of the tine, fatigue
python tillbridge/examples/chisel/make_figures.py out    # needs NumPy and Matplotlib
```

`out/Fig7_calibration.png` (and `.svg`) is Fig. 7 of the paper; `out/summary.txt` holds the
numbers. Finished steps are not repeated, so an interrupted run can be restarted with the same
command.

Seeds: the bed is generated with seed 11 (`BedOptions::seed`), the Latin-hypercube start design
and the batches of the calibration with seed 3 (`Calibration::set_seed`). All random numbers come
from the xoshiro256** generator, so the bed and the start design are identical on every platform.

The DEM runs of a granular flow are chaotic: a different compiler, compiler options or number of
threads changes the force history after a few tenths of a second, but not its mean and
statistics. The expected results below were obtained with GCC 13 (Linux, x86-64) and two threads;
with another platform the values must lie within the tolerances.

| Quantity (`out/summary.txt`) | Expected (GCC 13, 2 threads) | Tolerance on other platforms |
|---|---|---|
| reference draft, mean of 4 passes (`draft_reference`) | 238.495 N (s.d. 8.668 N) | ±5 % |
| reference vertical force (`vertical_reference`) | 25.489 N (s.d. 7.449 N) | ±10 N |
| estimate from both forces (`k_coh_estimate`, `mu_pw_estimate`) | 65.8 kPa, 0.573 | inside the 68 % contour of Fig. 7a |
| cohesion for the known friction μ_pw = 0.5 (`k_coh_calibrated`, `k_coh_low`, `k_coh_high`) | 85 kPa (75–95 kPa) | interval contains the true 90 kPa |
| number of calibration runs (`runs`) | 20 | exact |
| run time, two cores | 24–25 min (DEM step 4.1 ms) | depends on the processor |

The data of Fig. 7 are also given in `paper/figure_data/Fig7a_*.csv` and `Fig7b.csv`.

## Mouldboard plough (Section 4.1)

See `examples/plough/README.md` for the commands (DEM run, energy balance, wear, wrench, FE of
the shank in load cases A and D, fatigue, check against CalculiX). The DEM run takes about 35 min
on two cores.

Expected results of the repeated run with v2.1.0 (GCC 13, two threads) are listed in
`examples/plough/README.md`. The DEM run of the paper and the repeated run differ in their force
histories (chaotic flow); their mean draft differs by 5 %, and the fatigue lives of a single
block of 1.4 m differ markedly (Table 6 and Section 4.1 of the paper). The FE results of Table 5
(CalculiX check) use the mean loads given in `examples/plough/README.md` and are reproduced exactly.

## Data of the figures and tables

`paper/figure_data` holds the data of all figures and tables of the paper as CSV files (UTF-8,
comma separated, header in the first row). The same data are in the spreadsheet
`06_Figure_data_TillBridge.xlsx` submitted with the paper.

| File | Content |
|---|---|
| `Fig2a.csv` | oblique impact of a sphere on a wall: rebound velocities against incidence |
| `Fig2b.csv` | cantilever: tip deflection under mesh refinement, error to the extrapolated three-dimensional cantilever solution and difference to Timoshenko beam theory |
| `Fig2c.csv` | PCG iterations and times with AMG and Jacobi preconditioning |
| `Fig3_particles.csv` | positions and velocities of the particles of Fig. 3 |
| `Fig4a.csv` | energy balance of the plough body per bed and speed |
| `Fig4b.csv` | relative draft against ASABE D497 |
| `Fig5a.csv`, `Fig5b.csv` | contact pressure, frictional power density q_f and slip speed per facet at 1.0 and 2.5 m/s |
| `Fig6.csv` | wrench of load cases A and D and FE results |
| `Fig7a_runs.csv`, `Fig7a_contours.csv`, `Fig7a_map.csv`, `Fig7a_points.csv` | calibration runs, χ² contours, surrogate map, reference, estimate and calibrated cohesion |
| `Fig7b.csv` | forces of the calibration runs and of the four reference passes |
| `Table4.csv` | computing times |
| `Table5.csv` | verification cases |
| `Table5_shank_A.csv`, `Table5_shank_D.csv` | largest von Mises stress at the Gauss points of every element, TillBridge and CalculiX |
| `Table6.csv` | fatigue life for the detail categories (repeated run) |
| `Energy_closure.csv` | energy accounting of the repeated plough run and closure error against the time step |
