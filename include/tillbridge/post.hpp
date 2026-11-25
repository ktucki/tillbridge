// TillBridge - post-processing of DEM runs: energy balance, frictional power density (wear
// index), load transfer to the finite element model and fatigue.
// Copyright (c) 2026 Karol Tucki. MIT License.
#pragma once
#include <map>
#include <string>
#include <vector>

#include "tillbridge/dem.hpp"
#include "tillbridge/fem.hpp"
#include "tillbridge/tool.hpp"

namespace tb {

// Result files of a finished run (written by RunResult::write or by an earlier solver version).
struct RunFiles {
    std::vector<ForceSample> history;
    std::vector<TriangleLoad> triangles;   // Pt, Pn are zero if <prefix>_triW.csv is missing
    std::map<std::string, double> summary;
    Vec3 F_mean() const;
    static RunFiles read(const std::string& prefix);
};

// ---------------------------------------------------------------- energy balance and wear
// Work per metre of travel (J/m) in the averaging window. The draft work D equals the draft
// force; friction and normal damping at the tool come from the contact power; the rest is
// dissipated in the soil or stored as kinetic and potential energy of the slice
// (soil_and_slice = D - friction - normal).
struct EnergyBalance {
    double draft = 0;           // N = J/m
    double friction = 0;        // tangential contact power / travel speed
    double normal_damping = 0;  // normal contact power / travel speed
    double soil_and_slice = 0;
    double tangential_share = 0;   // share of the tangential contact forces in the draft force
};
EnergyBalance energy_balance(const RunFiles& r, const Tool& tool, double v_tool);

// Kinetic energy and lifting work of the soil leaving the tool, per metre of travel (J/m), from the
// departure records (<prefix>_depart.txt) of particles whose last tool contact ended with the tool
// between x0 and x1; lifting from the initial height of each body in the initial bed.
struct SliceEnergy {
    double kinetic = 0, lifting = 0;
    long particles = 0;
};
SliceEnergy slice_energy(const std::string& prefix, const Bed& initial, double x0, double x1, double g = 9.81);

// Per-triangle mean normal pressure p (Pa), frictional power density q_f (W/m^2) and mean slip
// speed v_s = P_t / |F_t| (m/s, for triangles with p > p_min). q_f is a relative wear index: by
// Archard's law the wear rate is K p v_s / H = K q_f / (mu H).
struct WearMap {
    std::vector<double> pressure, q_f, slip_speed;
    double total_frictional_power = 0;   // W
    double mean_pressure = 0;            // area-weighted over loaded triangles, Pa
    double median_slip_speed = 0;        // m/s, over loaded triangles
};
WearMap wear_map(const RunFiles& r, const Tool& tool, double p_min = 1000.0);
// Smooth a per-triangle field over neighbouring triangles (sharing a vertex), area weighted.
std::vector<double> smooth_on_tool(const Tool& tool, const std::vector<double>& f, int passes = 1);

// ---------------------------------------------------------------- load transfer
// Normal-force weighted centroid of the tool triangles (tool frame).
Vec3 centre_of_pressure(const Tool& tool, const std::vector<TriangleLoad>& loads);
// Force and moment about ref in the FE frame, with p_fe = R p_dem + t.
void wrench(const Vec3& F_dem, const Vec3& cop_dem, const double R[3][3], const Vec3& t, const Vec3& ref, Vec3& F, Vec3& M);
// Nodal forces on a node set that reproduce force F and moment M about ref (rigid spider, RBE3-like).
// Exact wrench of the facet forces of facets [first, end) about ref, in the frame of an FE model:
// p' = scale (R p + t), F' = R F; F = sum R F_k, M = sum (p'_k - ref) x R F_k (Eq. (8) of the paper).
void facet_wrench(const Tool& tool, const std::vector<TriangleLoad>& loads, size_t first, size_t end, const double R[3][3],
                  const Vec3& t, double scale, const Vec3& ref, Vec3& F, Vec3& M);
std::map<int64_t, Vec3> distribute(const Mesh& mesh, const std::vector<int64_t>& nodes, const Vec3& F, const Vec3& M, const Vec3& ref);
// |F(t)| with the tool between x0 and x1.
std::vector<double> force_magnitude(const std::vector<ForceSample>& h, double x0, double x1);

// ---------------------------------------------------------------- fatigue
struct Cycle {
    double range, mean, count;   // count 0.5 or 1
};
std::vector<double> turning_points(const std::vector<double>& x);
// Rainflow counting (ASTM E1049-85, three-point algorithm).
std::vector<Cycle> rainflow(const std::vector<double>& x);
// Endurance for a stress range (MPa), EN 1993-1-9 curve of detail category dsc, gamma_Mf = 1:
// slope 3 down to dsD = 0.737 dsc (5e6 cycles), slope 5 down to the cut-off dsL = 0.405 dsc (1e8).
double ec3_cycles(double ds, double dsc, bool cutoff = true);
double miner(const std::vector<Cycle>& c, double dsc, bool cutoff = true);
double equivalent_range(const std::vector<Cycle>& c, double m = 3.0);

struct FatigueResult {
    double cycles = 0, cycles_per_m = 0, max_range = 0, equivalent_range = 0, ratio_to_cutoff = 0;
    double damage_per_m = 0, life_km = 0, life_h = 0, damage_per_m_no_cutoff = 0, life_km_no_cutoff = 0;
};
// sigma(t) = sigma_ref |F(t)| / F_ref (linear elasticity), over a travel distance at speed v.
FatigueResult fatigue_life(const std::vector<double>& force_magnitude, double distance, double sigma_ref, double dsc, double speed,
                           double F_ref = 0);

}  // namespace tb
