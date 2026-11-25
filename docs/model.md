# Model

## Particles and integration

Bodies are rigid clusters of one or more spheres with a spherical moment of inertia I. Each
time step, contact forces and torques are summed per body, and translation and rotation are
advanced with the symplectic Euler scheme (velocity first, then position; the orientation
quaternion is updated with the exact rotation over dt). Bodies of material class 2 are fixed.
Neighbours are found on a uniform cell grid of size 2.05 r_max that is rebuilt every step, so
no neighbour list has to be kept up to date. The contact loop runs in parallel over spheres;
each contact is evaluated from both sides, so no atomic updates of the particles are needed and
results do not depend on the number of threads except through the order of the force sums on
the tool.

## Contact law

For two spheres with overlap δ, unit normal n, effective radius R* and effective mass m*,
effective modulus E* = [(1 − ν₁²)/E₁ + (1 − ν₂²)/E₂]⁻¹ and effective shear modulus
G* = [(2 − ν₁)/G₁ + (2 − ν₂)/G₂]⁻¹:

* normal force (Hertz with Tsuji damping):
  F_n = (4/3) E* (R*δ)^½ δ − 2 (5/6)^½ β (S_n m*)^½ v_n, with S_n = 2 E* (R*δ)^½,
  β = ln e / (ln² e + π²)^½, and F_n ≥ 0;
* cohesion (particle–particle only): F_c = k_coh π R* δ, subtracted from the normal force;
* tangential force (Mindlin with history): the tangential displacement u_t is accumulated from
  the relative velocity at the contact point and projected on the current tangent plane;
  F_t = −S_t u_t − 2 (5/6)^½ β (S_t m*)^½ v_t with S_t = 8 G* (R*δ)^½, limited by μ F_n
  (on sliding, u_t is reset to the value that gives the limit force);
* rolling resistance (particle–particle, type A elastic–plastic): a rolling moment from the
  relative rolling angular velocity with stiffness k_r = 2.25 k_n μ_r² R*², k_n = 2 E* (R*δ)^½,
  limited to μ_r R* F_n; twisting is not resisted.

The relative velocity is taken at the contact point, v = v_i + ω_i × (x_c − x_i) − (v_j +
ω_j × (x_c − x_j)), so rotation and sliding are coupled — the contact-law benchmarks in
`tests/test_dem.cpp` check this against rigid-body impact theory.

Walls (x = 0, x = L_x, y = 0, y = L_y, z = 0) and the tool are rigid and use the same law
with the wall or tool elastic constants, the particle–wall friction μ_pw, no cohesion and no
rolling resistance. A sphere touches the tool at the closest point of the nearest triangle;
its tangential history is kept when the contact point moves to a neighbouring triangle.

## Parallel bonds

With `bond_E` > 0, spheres of different bodies closer than (r_i + r_j)(1 + `bond_gap`) at t = 0
are joined by parallel bonds (Potyondy and Cundall 2004) of radius R_b = λ min(r_i, r_j),
normal stiffness k_n = E_b / (r_i + r_j) and shear stiffness k_n / k_r per unit area. Forces
and moments are accumulated incrementally; a bond breaks when the maximum tensile stress
−F_n/A + |M_b| R_b / I or the maximum shear stress |F_s|/A + |M_t| R_b / J exceeds its strength.
Bond damping is viscous with damping ratio `bond_zeta`.

## Recorded tool quantities

For each tool triangle the program averages over the averaging window the normal contact force,
the number of contacts, the force vector and the power dissipated at the contacts, split into
the tangential part P_t = −F_t · v_t (friction) and the normal part (damping). Their sums give
the energy balance of the draft work per metre of travel (`tillbridge energy`): D · 1 m =
friction + normal damping + work done in the soil. The frictional power density
q_f = P_t / A (W m⁻²) on each triangle is a relative wear index (`tillbridge wear`): by Archard's
law the wear rate is proportional to p v_s = q_f / μ_pw for sliding contacts. It is not calibrated
against measured material loss; a wear rate in m³ per metre of travel needs a wear coefficient
from tests.

## Energy accounting (version 2.1)

The work done in the soil can be obtained in two ways. `tillbridge energy` gives it as the
remainder of the draft work, D − (P_t + P_n)/v_tool. In addition, every run accounts directly for
the energy over the averaging window and writes it to `_summary.txt`:

* `E_tool` — work of the tool on the soil, the sum of F_c · v_tool over the tool contacts;
* `E_tool_contacts` — energy taken up at the tool contacts by the relative motion, the sum of
  −F_c · (v_i − v_tool) (the integral of P_t + P_n);
* `E_bulk` — energy taken up by all other contacts: particle–particle (including the rolling
  moments), walls and bonds, the sum of −(F_c · v_i + M_r · ω_i) taken from both bodies of each
  contact. It is the energy dissipated by friction and contact damping plus the change of the
  elastic and cohesive energy stored in the contacts;
* `dE_kinetic`, `dE_potential` — change of the kinetic and gravitational energy of all bodies
  between the start and the end of the window;
* `E_closure` = E_tool − E_tool_contacts − E_bulk − dE_kinetic − dE_potential.

The powers are evaluated with the same contact forces and contact-point velocities that drive
the bodies, so the bookkeeping is exact for each step; `E_closure` is the energy error of the
symplectic Euler integration. The solver has no numerical (global or local non-viscous) damping:
the only dissipation is the viscous contact damping of Eqs. (1)–(2), Coulomb friction and the
plastic part of the rolling resistance. `tillbridge energy` prints the directly accounted work
done in the soil, (E_bulk + dE_kinetic + dE_potential) per metre, beside the remainder, and the
closure error relative to the tool work.

## Stability

A run stops and is flagged as unstable when a particle exceeds `v_max` or, with `ke_growth`,
when a resting bed gains kinetic energy before the tool moves. The time step should be below
about 20 % of the Rayleigh time of the smallest particle.

## Finite elements

Linear elastic solids with 8-node and 20-node hexahedra (full 2×2×2 or 3×3×3 Gauss rules,
reduced 2×2×2 for imported SOLID186) and 10-node tetrahedra (4-point rule). The stiffness
matrix is assembled in compressed sparse row form and solved by the conjugate gradient method
preconditioned with smoothed-aggregation algebraic multigrid (geometric aggregates, the six
rigid-body modes as near-null space, damped-Jacobi prolongation smoothing, a dense Cholesky
factorisation on the coarse level) or with Jacobi.

## Load transfer

The coupling is one-way: the tool is rigid and moves along a prescribed straight path, its
deformation does not act back on the soil, and the loads are applied afterwards to a linear
elastic FE model. `tillbridge loads` reduces the window-mean force of every facet, acting at the
facet centroid, to the exact wrench about a reference point of the FE model (F = Σ R F_k,
M = Σ (p'_k − x_ref) × R F_k). The facet forces themselves are not mapped to FE nodes, because the
FE model is usually a different member (the shank) from the loaded surface (the body). The wrench
is applied as a distributed coupling (RBE3-like): nodal forces f_i = F/n + w × (x_i − c) on the
node set, with w = S⁻¹ M_c and S = Σ (|r_i|² I − r_i r_iᵀ), which reproduce F and the moment
exactly (to rounding). A rigid spider would add stiffness; this coupling does not.

## Fatigue

The tool force history in the averaging window gives a stress history σ(t) = σ_ref |F(t)| / F_ref
(linear elasticity, proportional loading: the direction and the point of application of the
load are taken as constant, so only its magnitude varies). σ_ref is the stress at the critical
location for the mean load, F_ref the mean of |F(t)| over the window. The window is one load
block; the damage of the block divided by the distance gives the damage per metre, and the life
in hours follows for continuous work at the travel speed of the run, L = 1 / (D_m/L_block · v ·
3600). The life therefore holds only for the operating spectrum represented by the block. Cycles are counted with the three-point rainflow
algorithm (ASTM E1049-85) and damage is summed with the Palmgren–Miner rule on the
EN 1993-1-9 S–N curve of the chosen detail category (slope 3 to 5·10⁶ cycles, slope 5 to the
cut-off at 10⁸ cycles, γ_Mf = 1); the life is also given without the cut-off.
