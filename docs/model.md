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
q_f = P_t / A on each triangle is a relative wear index (`tillbridge wear`).

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
factorisation on the coarse level) or with Jacobi. The soil load is applied as a distributed
coupling (RBE3-like): nodal forces on a node set that reproduce the force and the moment about a
reference point.

## Fatigue

The tool force history in the averaging window gives a stress history σ(t) = σ_ref |F(t)| / F_ref
(linear elasticity, proportional loading). Cycles are counted with the three-point rainflow
algorithm (ASTM E1049-85) and damage is summed with the Palmgren–Miner rule on the
EN 1993-1-9 S–N curve of the chosen detail category (slope 3 to 5·10⁶ cycles, slope 5 to the
cut-off at 10⁸ cycles, γ_Mf = 1); the life is also given without the cut-off.
