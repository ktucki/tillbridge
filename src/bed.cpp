// Copyright (c) 2026 Karol Tucki. MIT License.
#include "tillbridge/bed.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "tillbridge/textio.hpp"

namespace tb {

// ------------------------------------------------------------------ random numbers
static uint64_t splitmix(uint64_t& x) {
    uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
Rng::Rng(uint64_t seed) { for (auto& s : s_) s = splitmix(seed); }
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
uint64_t Rng::next() {
    uint64_t r = rotl(s_[1] * 5, 7) * 9, t = s_[1] << 17;
    s_[2] ^= s_[0]; s_[3] ^= s_[1]; s_[1] ^= s_[2]; s_[0] ^= s_[3]; s_[2] ^= t; s_[3] = rotl(s_[3], 45);
    return r;
}
double Rng::uniform(double a, double b) { return a + (b - a) * ((next() >> 11) * 0x1.0p-53); }
double Rng::normal() {
    if (has_spare_) { has_spare_ = false; return spare_; }
    double u, v, s;
    do { u = uniform(-1, 1); v = uniform(-1, 1); s = u * u + v * v; } while (s >= 1 || s == 0);
    double f = std::sqrt(-2 * std::log(s) / s);
    spare_ = v * f; has_spare_ = true;
    return u * f;
}
int Rng::integer(int lo, int hi) { return lo + int(next() % uint64_t(hi - lo)); }

// ------------------------------------------------------------------ bed
double Body::top() const {
    double t = -1e300;
    for (auto& s : spheres) t = std::max(t, x.z + s.offset.z + s.r);
    return t;
}
size_t Bed::n_spheres() const { size_t n = 0; for (auto& b : bodies) n += b.spheres.size(); return n; }
double Bed::mass() const { double m = 0; for (auto& b : bodies) m += b.m; return m; }

Bed Bed::read(const std::string& path) {
    TextReader r(path);
    std::vector<double> h = r.numbers();
    if (h.size() < 2) r.fail("bad bed header");
    long NB = long(h[0]), NS = long(h[1]);
    if (NB <= 0 || NS < NB) r.fail("bad bed header");
    Bed bed; bed.bodies.resize(NB);
    long k = 0;
    for (auto& b : bed.bodies) {
        auto v = r.numbers();
        size_t n = v.size();
        if (n < 6 || (n > 7 && n != 10 && n != 13)) r.fail("bad body record");
        b.x = {v[0], v[1], v[2]}; b.m = v[3]; b.I = v[4];
        int ns = int(v[5]);
        b.mat = n >= 7 ? int(v[6]) : 0;
        if (n >= 10) b.v = {v[7], v[8], v[9]};
        if (n >= 13) b.w = {v[10], v[11], v[12]};
        if (ns < 1 || b.m <= 0 || b.I <= 0 || k + ns > NS) r.fail("bad body record");
        b.spheres.resize(ns);
        for (auto& s : b.spheres) {
            auto q = r.numbers();
            if (q.size() < 4 || q[3] <= 0) r.fail("bad sphere record");
            s.offset = {q[0], q[1], q[2]}; s.r = q[3];
        }
        k += ns;
    }
    if (k != NS) r.fail("sphere count differs from header");
    return bed;
}

void Bed::write(const std::string& path, bool with_velocity) const {
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) throw std::runtime_error("cannot write " + path);
    std::fprintf(f, "%zu %zu\n", size(), n_spheres());
    for (auto& b : bodies) {
        std::fprintf(f, "%.6f %.6f %.6f %.6e %.6e %zu %d", b.x.x, b.x.y, b.x.z, b.m, b.I, b.spheres.size(), b.mat);
        if (with_velocity)
            std::fprintf(f, " %.9e %.9e %.9e %.9e %.9e %.9e", b.v.x, b.v.y, b.v.z, b.w.x, b.w.y, b.w.z);
        std::fprintf(f, "\n");
        for (auto& s : b.spheres) std::fprintf(f, "%.6f %.6f %.6f %.6f\n", s.offset.x, s.offset.y, s.offset.z, s.r);
    }
    std::fclose(f);
}

Bed Bed::select(const std::function<bool(const Body&)>& keep) const {
    Bed out;
    for (auto& b : bodies) if (keep(b)) out.bodies.push_back(b);
    return out;
}

// ------------------------------------------------------------------ generation
static double union_volume(const std::vector<Sphere>& s, double ext, const std::vector<Vec3>& mc) {
    long in = 0;
    for (auto& p0 : mc) {
        Vec3 p = p0 * ext;
        for (auto& q : s) { Vec3 d = p - q.offset; if (dot(d, d) <= q.r * q.r) { in++; break; } }
    }
    return double(in) / mc.size() * std::pow(2 * ext, 3);
}

std::vector<Sphere> clod_shape(double R, Rng& rng, double* bounding, int n_mc, double max_aspect) {
    std::vector<Vec3> mc(n_mc);
    for (auto& p : mc) p = {rng.uniform(-1, 1), rng.uniform(-1, 1), rng.uniform(-1, 1)};
    for (;;) {
        int k = rng.integer(3, 5);
        std::vector<Sphere> s(k);
        Vec3 mean;
        for (auto& q : s) {
            q.r = R * rng.uniform(0.55, 0.75);
            Vec3 c{rng.normal(), rng.normal(), rng.normal()};
            c = c * (1.0 / norm(c));
            q.offset = c * ((R - q.r) * rng.uniform(0.7, 1.0));
            mean += q.offset;
        }
        mean = mean * (1.0 / k);
        double ext = 0;
        for (auto& q : s) { q.offset -= mean; ext = std::max(ext, norm(q.offset) + q.r); }
        double sc = R / ext;
        for (auto& q : s) { q.offset *= sc; q.r *= sc; }
        double V = union_volume(s, R, mc);
        double f = std::cbrt(4.0 / 3.0 * kPi * R * R * R / V);
        ext = 0;
        for (auto& q : s) { q.offset *= f; q.r *= f; ext = std::max(ext, norm(q.offset) + q.r); }
        if (ext / R <= max_aspect) { if (bounding) *bounding = ext; return s; }
    }
}

static Vec3 rotate_random(const Vec3& p, const double Rm[3][3]) {
    return {Rm[0][0] * p.x + Rm[0][1] * p.y + Rm[0][2] * p.z, Rm[1][0] * p.x + Rm[1][1] * p.y + Rm[1][2] * p.z,
            Rm[2][0] * p.x + Rm[2][1] * p.y + Rm[2][2] * p.z};
}

Bed generate_bed(int n, const BedOptions& o) {
    if (n <= 0) throw std::invalid_argument("number of particles must be positive");
    Rng rng(o.seed);
    struct Part { double ext, m; std::vector<Sphere> s; };
    std::vector<Part> parts;
    for (int i = 0; i < n; i++) {
        double R = o.radius * rng.uniform(1 - o.spread, 1 + o.spread);
        double m = o.density * 4.0 / 3.0 * kPi * R * R * R;
        if (!o.clods) { parts.push_back({R, m, {Sphere{Vec3{}, R}}}); continue; }
        double ext;
        auto s = clod_shape(R, rng, &ext);
        Quat q{rng.normal(), rng.normal(), rng.normal(), rng.normal()};
        normalize(q);
        double w = q.w, a = q.x, b = q.y, d = q.z;
        double Rm[3][3] = {{1 - 2 * (b * b + d * d), 2 * (a * b - d * w), 2 * (a * d + b * w)},
                           {2 * (a * b + d * w), 1 - 2 * (a * a + d * d), 2 * (b * d - a * w)},
                           {2 * (a * d - b * w), 2 * (b * d + a * w), 1 - 2 * (a * a + b * b)}};
        for (auto& sp : s) sp.offset = rotate_random(sp.offset, Rm);
        parts.push_back({ext, m, s});
    }
    double emax = 0; for (auto& p : parts) emax = std::max(emax, p.ext);
    double sp = 2 * emax + 0.001;
    std::vector<double> xs, ys;
    for (double x = sp / 2 + 0.005; x < o.Lx - sp / 2; x += sp) xs.push_back(x);
    for (double y = sp / 2 + 0.005; y < o.Ly - sp / 2; y += sp) ys.push_back(y);
    if (xs.empty() || ys.empty()) throw std::invalid_argument("box too small for the particles");
    std::vector<Vec3> cols; for (double y : ys) for (double x : xs) cols.push_back({x, y, 0});
    std::vector<double> ztop(cols.size(), 0.0);
    Bed bed; bed.bodies.resize(n);
    for (int i = 0; i < n; i++) {
        auto& P = parts[i]; size_t j = i % cols.size();
        double z = ztop[j] + P.ext + 0.001; ztop[j] = z + P.ext;
        Body& b = bed.bodies[i];
        b.x = {cols[j].x + rng.uniform(-1e-3, 1e-3), cols[j].y + rng.uniform(-1e-3, 1e-3), z};
        double Req = std::cbrt(3 * P.m / (4 * kPi * o.density));
        b.m = P.m; b.I = 0.4 * P.m * Req * Req; b.spheres = P.s;
    }
    return bed;
}

int particles_for_depth(double Lx, double Ly, double depth, double radius, double density, double bulk) {
    double m = density * 4.0 / 3.0 * kPi * radius * radius * radius;
    return int(std::ceil(bulk * Lx * Ly * depth / m * 1.08));
}

Bed trim(const Bed& b, double depth, double tol) {
    return b.select([&](const Body& x) { return x.top() <= depth + tol; });
}

Bed cut_furrow(const Bed& b, double y_edge, double z_bottom, double angle_deg, bool low) {
    double t = 1.0 / std::tan(angle_deg * kPi / 180.0);
    return b.select([&](const Body& x) {
        double y = x.x.y, z = x.x.z;
        bool cut = low ? (z > z_bottom && y < y_edge - (z - z_bottom) * t) : (z > z_bottom && y > y_edge + (z - z_bottom) * t);
        return !cut;
    });
}

void fix_below(Bed& b, double z) { for (auto& x : b.bodies) if (x.x.z < z) x.mat = kFixed; }

std::map<std::string, double> soil_fate(const Bed& a, const Bed& f, const std::function<bool(const Vec3&)>& region,
                                        const std::vector<std::pair<std::string, std::function<bool(const Vec3&)>>>& zones) {
    if (a.size() != f.size()) throw std::invalid_argument("beds differ in size");
    std::map<std::string, double> out;
    std::vector<char> used(a.size(), 0);
    double M = 0;
    for (size_t i = 0; i < a.size(); i++) if (region(a.bodies[i].x)) M += a.bodies[i].m;
    if (M <= 0) throw std::invalid_argument("region is empty");
    double rest = 0;
    for (auto& z : zones) {
        double m = 0;
        for (size_t i = 0; i < a.size(); i++)
            if (!used[i] && region(a.bodies[i].x) && z.second(f.bodies[i].x)) { used[i] = 1; m += a.bodies[i].m; }
        out[z.first] = m / M;
    }
    for (size_t i = 0; i < a.size(); i++) if (!used[i] && region(a.bodies[i].x)) rest += a.bodies[i].m;
    out["other"] = rest / M;
    return out;
}

}  // namespace tb
