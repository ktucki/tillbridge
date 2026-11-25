# Mouldboard plough: from the soil bin to the fatigue life of the shank

The example of Section 4.1 of the paper. All inputs are in this folder:

| File | Content |
|---|---|
| `bed.txt` | settled soil bin 5.0 × 1.5 m, 47 928 spheres |
| `tool.txt` | share, mouldboard and frog, 7184 triangles (m); facets 0–7143 belong to the body |
| `params.txt` | calibrated DEM parameters, travel speed 1.5 m/s |
| `shank.msh` | shank (mm), 166 095 nodes, 20-node hexahedra and 10-node tetrahedra, Gmsh 2.2; node groups `CLAMP` (frame flange) and `LOADFACE` (frog flange) |
| `compare_calculix.py` | comparison of the FE solution with CalculiX on the same mesh and nodal loads |

The DEM frame (m) is turned into the frame of the shank model (mm) by p' = 1000 (R p + T) with

```
R = [ 0.5090414158  0.8607420270  0
      0             0            -1
     -0.8607420270  0.5090414158  0 ]
T = [-0.6591927285  0.6974405929 -1.4319915082]
```

## Commands

```
# 1. DEM run (about 30 min on two cores)
tillbridge run bed.txt tool.txt P15 --params params.txt --threads 2

# 2. energy balance and wear index over the steady part of the run
tillbridge energy P15 tool.txt --speed 1.5 --bed bed.txt --x0 2.02 --x1 3.42
tillbridge wear P15 tool.txt --smooth 1 --out P15_wear.csv

# 3. exact wrench of the body facets in the frame of the shank, about the reference point (mm)
tillbridge loads P15 tool.txt --facets 0 7144 --scale 1000 \
    --rotation 0.5090414158 0.8607420270 0 0 0 -1 -0.8607420270 0.5090414158 0 \
    --translation -0.6591927285 0.6974405929 -1.4319915082 \
    --ref 109.65347574 516.80816546 -406.45131372

# 4. load case A: the exact wrench (force F and moment M_A from step 3)
tillbridge fem shank.msh --E 200000 --nu 0.3 --fixed CLAMP --loaded LOADFACE \
    --force -1304.41 347.594 3957.55 --moment -254573 188732 -63871.4 \
    --ref 109.65347574 516.80816546 -406.45131372 --vtk caseA.vtk --ccx caseA.inp

# 5. load case D: the same force at one third of the working depth (free couple kept)
tillbridge fem shank.msh --E 200000 --nu 0.3 --fixed CLAMP --loaded LOADFACE \
    --force -1304.41 347.594 3957.55 --moment -587.625 188732 19842.3 \
    --ref 109.65347574 516.80816546 -406.45131372 --vtk caseD.vtk --ccx caseD.inp

# 6. fatigue life with the largest von Mises stress in the shank tube of each case
tillbridge fatigue P15 --x0 2.02 --x1 3.42 --speed 1.5 --detail 71 --sigma-ref 37.1
tillbridge fatigue P15 --x0 2.02 --x1 3.42 --speed 1.5 --detail 71 --sigma-ref 65.5

# 7. optional check with CalculiX (ccx 2.22): ccx -i caseA, then
python3 compare_calculix.py caseA caseA.vtk
```

## Results (paper, Section 4.1 and Table 5)

| | case A | case D |
|---|---|---|
| largest von Mises stress in the tube (MPa) | 37.1 | 65.5 |
| largest stress range (MPa) | 30.2 | 53.4 |
| fatigue life (h of ploughing), detail category 71 | about 40 000 | about 1400 |
| displacements against CalculiX (relative norm) | 1.6 × 10⁻⁶ | 1.2 × 10⁻⁶ |
| Gauss-point stresses against CalculiX (relative norm) | 1.2 × 10⁻⁷ | 1.2 × 10⁻⁷ |

The FE cases use the mean loads of the run of the paper, given above, so that steps 4–7
reproduce Table 5 exactly. The DEM run is chaotic: a different compiler or number of threads
changes the force history after a few tenths of a second, but not its mean and statistics.

## Repeated run with version 2.1.0 (GCC 13, two threads)

| | value |
|---|---|
| draft (mean, steady part) | 5227 N (paper: 4960 and 5393 N on two beds) |
| wrench about the reference point, step 3 | F = (−1519.0, 368.7, 4059.6) N, M = (−248 620, 143 750, −90 821) N mm |
| sum of facet forces against the mean resultant | relative difference 2 × 10⁻⁸ |
| work done in the soil: remainder / direct accounting | 3544 / 3742 J/m |
| energy closure error | −198 J/m, −3.8 % of the draft work (first order in the time step) |
| fatigue block | 1.4 m, 0.93 s, force every 0.01 s, 27 cycles (19.3 per metre) |
| case D moment (same force at one third of the depth) | M = (2964.8, 143 750, 3313.8) N mm |
| von Mises stress at the critical node, case A / case D | 39.5 / 68.1 MPa (paper: 37.1 / 65.5 MPa) |
| largest range, case A / case D | 20.8 / 35.8 MPa |
| life, detail category 71, case A / case D | unlimited (below the cut-off) / about 5300 h |
| critical node of the tube (both cases) | node 21456 of shank.msh, (95.7, 337.0, −219.9) mm |
| reactions at CLAMP against the applied force and moment | equal to the six printed digits |

The lives depend on the few largest cycles of a short block: the run of the paper gave about
40 000 h (case A) and 1400 h (case D); see Section 4.1 and Table 6 of the paper.
