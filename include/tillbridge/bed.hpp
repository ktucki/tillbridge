// TillBridge - soil beds of rigid multi-sphere bodies.
// Copyright (c) 2026 Karol Tucki. MIT License.
#pragma once
#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "tillbridge/vec.hpp"

namespace tb {

// One sphere of a body: centre relative to the body centre (box axes) and radius.
struct Sphere {
    Vec3 offset;
    double r = 0;
};

// Material classes of bodies.
enum Material : int { kSoil = 0, kStone = 1, kFixed = 2 };

struct Body {
    Vec3 x;                  // centre
    double m = 0;            // mass
    double I = 0;            // moment of inertia (spherical)
    std::vector<Sphere> spheres;
    int mat = kSoil;
    Vec3 v, w;               // initial velocity and angular velocity
    double top() const;      // highest point
};

// A soil bed. File format (text, SI units):
//   NB NS
//   x y z m I ns [mat [vx vy vz [wx wy wz]]]
//   dx dy dz r      (ns lines)
struct Bed {
    std::vector<Body> bodies;

    size_t size() const { return bodies.size(); }
    size_t n_spheres() const;
    double mass() const;

    static Bed read(const std::string& path);
    // Velocities are written only when with_velocity is true.
    void write(const std::string& path, bool with_velocity = false) const;
    Bed select(const std::function<bool(const Body&)>& keep) const;
};

struct BedOptions {
    double Lx = 1.0, Ly = 1.0;
    double radius = 0.015;      // volume-equivalent radius
    double spread = 0.15;       // uniform variation of the radius, +-
    bool clods = false;         // false: spheres; true: 3-4 overlapping spheres
    double density = 2600.0;    // particle (solid) density
    uint64_t seed = 1;
};

// Particles stacked in columns above the floor, to be settled under gravity.
Bed generate_bed(int n, const BedOptions& o);
// Number of particles that settle to about `depth` at the given bulk density.
int particles_for_depth(double Lx, double Ly, double depth, double radius, double density = 2600.0,
                        double bulk_density = 1500.0);
// Keep bodies whose highest point is below depth + tol.
Bed trim(const Bed& b, double depth, double tol = 0.005);
// Remove the soil of the previous furrow: bodies above z_bottom on the low-y side of a wall
// that starts at y_edge on the furrow bottom and rises at angle_deg.
Bed cut_furrow(const Bed& b, double y_edge, double z_bottom, double angle_deg = 45.0, bool low_side = true);
// Mark bodies whose centre is below z as fixed (material class 2, anchored subsoil).
void fix_below(Bed& b, double z);

// Irregular clod of 3-4 spheres whose union has the volume of a sphere of radius R
// (union volume by Monte Carlo). Returns the spheres and the bounding radius.
struct Rng;
std::vector<Sphere> clod_shape(double R, Rng& rng, double* bounding_radius, int n_mc = 20000, double max_aspect = 1.30);

// Mass fractions of the soil that started in `region` and ended in each zone (bodies matched by index).
std::map<std::string, double> soil_fate(const Bed& initial, const Bed& final_state,
                                        const std::function<bool(const Vec3&)>& region,
                                        const std::vector<std::pair<std::string, std::function<bool(const Vec3&)>>>& zones);

// Deterministic random numbers (xoshiro256**), identical on all platforms.
struct Rng {
    explicit Rng(uint64_t seed);
    uint64_t next();
    double uniform(double a = 0.0, double b = 1.0);
    double normal();
    int integer(int lo, int hi_exclusive);
  private:
    uint64_t s_[4];
    bool has_spare_ = false;
    double spare_ = 0;
};

}  // namespace tb
