// TillBridge - discrete element solver for soil-tool interaction.
// Copyright (c) 2026 Karol Tucki. MIT License.
//
// Rigid multi-sphere bodies in a box with a rigid triangulated tool moving at constant
// velocity along +x. Contact law: Hertz-Mindlin with viscous damping from the coefficient of
// restitution (Tsuji et al. 1992), Coulomb friction with tangential history, linear cohesion
// k_coh * pi * R* * overlap, and rolling resistance (elastic-plastic spring, model C of
// Ai et al. 2011 without the viscous term). Optional parallel bonds (Potyondy and Cundall 2004).
// Velocities are evaluated at the contact point, so rotation and sliding are coupled.
// The power dissipated at every tool contact is accumulated per triangle (tangential and
// normal part), which gives the energy balance and the frictional power density (wear index).
#pragma once
#include <map>
#include <string>
#include <vector>

#include "tillbridge/bed.hpp"
#include "tillbridge/params.hpp"
#include "tillbridge/tool.hpp"

namespace tb {

constexpr const char* kVersion = "2.1.0";

struct ForceSample {
    double t, Fx, Fy, Fz, Ekin, x_tool;
};

// Per-triangle means over the averaging window.
struct TriangleLoad {
    double Fn = 0;          // normal contact force (compression positive), N
    double contacts = 0;    // number of contacts
    Vec3 F;                 // force of the soil on the triangle, N
    double Pt = 0;          // power dissipated by tangential forces (friction), W
    double Pn = 0;          // power dissipated by normal forces (damping), W
};

// Particle leaving the tool: position, velocity and tool position at the end of its last contact.
struct Departure {
    long sphere, body;
    Vec3 x, v;
    double x_tool;
};

struct Bond {
    long body_i, body_j;
    bool alive;
};

// Energy accounting over the averaging window (J): work of the tool on the soil, energy taken up
// at the tool contacts by the relative motion (friction and normal damping, Eq. (5) of the paper),
// energy taken up by all other contacts (particle-particle, walls, rolling resistance, bonds:
// dissipation plus the change of the energy stored in the contacts), and the change of the
// kinetic and potential energy of the soil. closure() is the part of the tool work that is not
// accounted for, i.e. the error of the time integration.
struct EnergyAccount {
    double tool = 0, tool_contacts = 0, bulk = 0, kinetic = 0, potential = 0, distance = 0;
    double closure() const { return tool - tool_contacts - bulk - kinetic - potential; }
};

struct RunResult {
    bool unstable = false;
    double t_unstable = 0, v_max_found = 0, ke_unstable = 0;
    std::vector<ForceSample> history;
    std::vector<TriangleLoad> triangles;
    std::vector<Departure> departures;
    std::vector<Bond> bonds;
    Bed state;                       // final positions (sphere offsets rotated) and velocities
    long steps = 0, n_avg = 0, bonds_broken_tension = 0, bonds_broken_shear = 0;
    double cpu_s = 0, mean_contacts_per_step = 0;
    Vec3 F_mean;                     // mean force of the soil on the tool in the averaging window
    EnergyAccount energy;            // energy accounting over the averaging window

    double draft() const { return -F_mean.x; }
    // Writes <prefix>_force.csv, _tri.csv, _triW.csv, _state.txt, _vel.txt, _depart.txt,
    // [_bonds.txt] and _summary.txt (the formats of docs/file_formats.md).
    void write(const std::string& prefix) const;
    std::map<std::string, double> summary() const;
};

class Simulation {
  public:
    Simulation(const Params& p, const Bed& bed, const Tool& tool);
    ~Simulation();
    Simulation(const Simulation&) = delete;
    Simulation& operator=(const Simulation&) = delete;
    // Runs to t_end (or until the run is found unstable).
    RunResult run(bool verbose = false);

  private:
    struct Impl;
    Impl* d_;
};

// Convenience: run and write the outputs under p.out. Returns the result.
RunResult run_dem(const Params& p, const Bed& bed, const Tool& tool, bool write = true, bool verbose = false);

// Contact force of the Hertz-Mindlin law (exposed for testing). n points from the partner to the
// particle; vrel is the particle velocity relative to the partner at the contact point; ut is the
// tangential spring history (updated). Returns the force on the particle.
Vec3 hertz_mindlin(double overlap, const Vec3& n, const Vec3& vrel, double Reff, double meff, double Eeff,
                   double Geff, double mu, double beta, double dt, double k_coh, double area, Vec3& ut);
double damping_beta(double restitution);
double effective_modulus(double E1, double nu1, double E2, double nu2);
double effective_shear_modulus(double E1, double nu1, double E2, double nu2);
// Closest point of triangle (a, b, c) to p (Ericson, Real-Time Collision Detection 5.1.5).
Vec3 closest_on_triangle(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c);
// Rayleigh time of a sphere of radius R (s).
double rayleigh_time(double R, double density, double E, double nu);

}  // namespace tb
