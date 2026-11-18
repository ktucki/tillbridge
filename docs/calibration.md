# Calibration

`tillbridge calibrate CONFIG` searches DEM parameters with Gaussian-process surrogates
(Matérn 5/2 kernel with one length scale per parameter, hyperparameters by maximum marginal
likelihood) and expected improvement. Every simulation is expensive, so the search starts from
a Latin hypercube design and then proposes a few runs at a time.

## Two problem types

* **Target matching** (one or more `target` lines): minimise
  Σ w ((y − y*) / y*)² over the measured quantities y* (for example the draft force and the
  vertical force measured in a soil bin).
* **Constrained** (`objective` and `constraint` lines): minimise one quantity while others stay
  in given intervals; the acquisition is the expected improvement times the probability that
  all constraints are met (Gardner et al. 2014; Gelbart et al. 2014).

## Configuration file

One item per line, `#` starts a comment.

| Item | Meaning |
|---|---|
| `bed FILE` | bed file (settled and prepared) |
| `tool FILE` | tool triangle file |
| `params FILE` | parameter file with the fixed parameters |
| `set KEY VALUE` | change one fixed parameter (may be repeated) |
| `param KEY LOW HIGH` | a parameter to calibrate and its bounds (may be repeated) |
| `target QUANTITY VALUE [WEIGHT]` | a measured quantity (target matching) |
| `objective QUANTITY` | the quantity to minimise (constrained problem) |
| `constraint QUANTITY LOW HIGH` | an interval for a quantity (constrained problem) |
| `init N` | size of the start design (default 8) |
| `iter N` | number of Bayesian steps (default 10) |
| `batch N` | runs proposed per step (kriging believer, default 1) |
| `seed N` | random seed of the start design |
| `log FILE` | evaluation log (default `calibration.log`) |
| `out PREFIX` | prefix of the run files (`PREFIX_cal_000`, ...) |

The quantities are the keys of the run summary — `Fx_mean`, `Fy_mean`, `Fz_mean`,
`mean_contacts_per_step`, `bonds_broken_tension`, ... — and `draft` (= −`Fx_mean`). A run that
becomes numerically unstable is scored as zero load, so the search moves away from it.

Example (the draft and the vertical force of a chisel tine; the cohesion and the soil–steel
friction are calibrated, the particle friction is fixed — the draft alone cannot separate it
from the cohesion):

```
bed     bed.txt
tool    tine.txt
params  base.txt
set     t_end 0.63
set     avg_x0 0.30
set     avg_x1 0.75
set     mu_pp 0.55
param   k_coh 0 1.8e5
param   mu_pw 0.2 0.9
target  draft 238.0
target  Fz_mean -25.0
init    8
iter    6
batch   2
```

## Scattered results

A DEM run is not exactly repeatable in the sense that a small change of a parameter or of the
tool position rearranges the particles: in the chisel example the draft of single runs with the
same parameters scatters by about 3 % and the vertical force by about 20 %. The best single run
can therefore be a lucky one. The program reports two answers:

* the estimate (the parameter names in the output, with `predicted_` quantities): the minimum
  of the misfit, or of the objective under the constraints, computed from Gaussian-process
  means of the individual quantities; the processes include a noise term, so they average
  neighbouring runs. This is the recommended result;
* the best single run (`best_run_` entries).

Run the solver once more with the estimated parameters to obtain their result files. Measured
targets should be means of replicate passes, as in a soil bin.

## Log and restart

Each evaluation is appended to the log as one line

```
label source x_1 ... x_d | key=value key=value ...
```

where `source` is `lhs` for the start design and `boN` for Bayesian step N. When the program
is started again with the same log, the evaluations already in it are read, the start design is
skipped and the search continues. Results of runs made elsewhere (for example with a different
solver or by hand) can be added as lines with any label and source.

## Library

The same procedure is available as `tb::Calibration` (`include/tillbridge/calibrate.hpp`) with
any user function that maps a parameter vector to named results, so other quantities — soil
throw, furrow profile, forces from several speeds — can be matched as well. See
`examples/chisel/chisel_example.cpp`.
