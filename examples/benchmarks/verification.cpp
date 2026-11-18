// Verification table: every test case with its reference solution and the largest relative
// difference. Writes verification.csv (case, reference, quantity, TillBridge, reference value,
// relative difference).
// Copyright (c) 2026 Karol Tucki. MIT License.
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "tillbridge/dem.hpp"
#include "tillbridge/fem.hpp"
#include "tillbridge/post.hpp"

using namespace tb;

namespace {
const double R = 0.015, RHO = 2600.0, M = RHO * 4.0 / 3.0 * kPi * R * R * R, I0 = 0.4 * M * R * R;
std::FILE* out;
void row(const char* c, const char* ref, const char* q, double sim, double exact) {
    double d = exact != 0 ? sim / exact - 1 : sim;
    std::fprintf(out, "\"%s\",\"%s\",\"%s\",%.6g,%.6g,%.3e\n", c, ref, q, sim, exact, d);
    std::printf("%-44s %-28s sim %-12.6g ref %-12.6g diff %+.2e\n", c, q, sim, exact, d);
}
Bed spheres(const std::vector<Vec3>& x, const std::vector<Vec3>& v = {}, const std::vector<Vec3>& w = {}) {
    Bed b;
    for (size_t i = 0; i < x.size(); i++) {
        Body B;
        B.x = x[i]; B.m = M; B.I = I0; B.spheres = {Sphere{Vec3{}, R}};
        if (i < v.size()) B.v = v[i];
        if (i < w.size()) B.w = w[i];
        b.bodies.push_back(B);
    }
    return b;
}
Params base() {
    Params p;
    p.E = 2.6e6; p.nu = 0.3; p.mu_pp = 0; p.mu_pw = 0; p.e = 0.5; p.mu_r = 0; p.k_coh = 0; p.dt = 1e-6; p.t_tool = 99; p.v_tool = 0;
    p.Lx = 2.0; p.Ly = 0.4; p.Lz = 0.4; p.out_dt = 1e-4;
    return p;
}
}  // namespace

int main() {
    out = std::fopen("verification.csv", "w");
    std::fprintf(out, "case,reference,quantity,tillbridge,reference_value,relative_difference\n");
    {   // static overlap
        Params p = base(); p.e = 0.1; p.t_end = 0.3; p.dt = 2e-6;
        RunResult r = run_dem(p, spheres({{0.5, 0.2, R}}), {}, false);
        double Ee = effective_modulus(p.E, p.nu, p.E_wall, p.nu_wall);
        row("sphere resting on the floor", "Hertz", "overlap", R - r.state.bodies[0].x.z, std::pow(M * p.g / (4.0 / 3.0 * Ee * std::sqrt(R)), 2.0 / 3.0));
    }
    {   // restitution independent of speed
        double e1 = 0, e2 = 0;
        for (double vn : {1.5, 3.0}) {
            Params p = base(); p.g = 0; p.e = 0.8; p.t_end = 0.02;
            RunResult r = run_dem(p, spheres({{0.5, 0.2, R + 0.001}}, {{0, 0, -vn}}), {}, false);
            (vn < 2 ? e1 : e2) = r.state.bodies[0].v.z / vn;
        }
        row("normal impact, 1.5 and 3 m/s", "Tsuji et al. damping", "rebound ratio at 3 vs 1.5 m/s", e2, e1);
    }
    {   // oblique impact, sliding
        double mu = 0.3, vn = 1.0, vt = 3.0;
        Params p = base(); p.E = 2.6e8; p.mu_pw = mu; p.g = 0; p.t_end = 0.02;
        RunResult r = run_dem(p, spheres({{0.5, 0.2, R + 0.002}}, {{vt, 0, -vn}}), {}, false);
        const Body& b = r.state.bodies[0];
        double e = b.v.z / vn;
        row("oblique impact on a wall, sliding", "Chung and Ooi, test 6", "tangential velocity", b.v.x, vt - mu * (1 + e) * vn);
        row("oblique impact on a wall, sliding", "Chung and Ooi, test 6", "angular velocity", b.w.y, 2.5 * mu * (1 + e) * vn / R);
    }
    {   // spinning sphere set down
        double w0 = 100.0;
        Params p = base(); p.mu_pw = 0.5; p.e = 0.3; p.t_end = 0.4; p.dt = 2e-6;
        RunResult r = run_dem(p, spheres({{0.5, 0.2, R - 1e-5}}, {}, {{0, w0, 0}}), {}, false);
        row("spinning sphere set down on a floor", "rigid-body dynamics", "final velocity", r.state.bodies[0].v.x, 2.0 / 7.0 * R * w0);
    }
    {   // tool plate drags a sphere
        double vtool = 0.5, z = 0.2;
        Tool t;
        Vec3 P0{-1.0, 0.0, z}, P1{1.5, 0.0, z}, P2{1.5, 0.4, z}, P3{-1.0, 0.4, z};
        t.tri = {{P0, P1, P2}, {P0, P2, P3}};
        t = t.refined(0.05);
        Params p = base(); p.mu_pw = 0.5; p.e = 0.3; p.t_tool = 0; p.v_tool = vtool; p.t_end = 0.3; p.dt = 2e-6; p.avg_x0 = -1; p.avg_x1 = 2;
        RunResult r = run_dem(p, spheres({{0.5, 0.2, z + R - 4e-5}}), t, false);
        row("sphere on a moving tool plate", "rigid-body dynamics", "final velocity", r.state.bodies[0].v.x, 2.0 / 7.0 * vtool);
    }
    {   // binary oblique collision
        double mu = 0.2, vn = 1.0, vt = 1.6;
        Params p = base(); p.E = 2.6e11; p.mu_pp = mu; p.e = 0.6; p.g = 0; p.t_end = 0.0008; p.dt = 1e-8;
        RunResult r = run_dem(p, spheres({{0.5, 0.2, 0.2}, {0.5 + 2 * R + 1e-7, 0.2, 0.2}}, {{vn / 2, vt / 2, 0}, {-vn / 2, -vt / 2, 0}}), {}, false);
        const Body &a = r.state.bodies[0], &b = r.state.bodies[1];
        double Jn = M * (b.v.x - (-vn / 2));
        row("oblique collision of two spheres", "Chung and Ooi, test 7", "angular velocity", std::fabs(a.w.z), 2.5 * mu * Jn / (M * R));
    }
    {   // FE: cantilever, hex20
        Vec3 L{400.0, 20.0, 40.0};
        Mesh m = box_mesh(L, 60, 6, 12, {}, 2);
        std::vector<int64_t> fixed, end;
        for (size_t i = 0; i < m.n_nodes(); i++) {
            if (m.X[i].x < 1e-9) fixed.push_back(int64_t(i));
            if (m.X[i].x > L.x - 1e-9) end.push_back(int64_t(i));
        }
        std::map<int64_t, Vec3> loads;
        for (auto i : end) loads[i] = {0, 0, -1000.0 / end.size()};
        SolverOptions o; o.rtol = 1e-10;
        Solution s = solve_static(m, 2.1e5, 0.3, fixed, loads, o);
        double w = 0;
        for (auto i : end) w -= s.U[i].z / end.size();
        double E = 2.1e5, b = 20, h = 40, Ll = 400, P = 1000, I = b * h * h * h / 12, G = E / 2.6, A = b * h;
        row("cantilever, 20-node hexahedra", "Timoshenko beam", "tip deflection", w, P * Ll * Ll * Ll / (3 * E * I) + P * Ll / (5.0 / 6.0 * G * A));
        row("cantilever, 20-node hexahedra", "equilibrium", "reaction", s.reaction.z, P);
    }
    {   // fatigue
        std::map<double, double> count;
        for (const Cycle& c : rainflow({-2, 1, -3, 5, -1, 3, -4, 4, -2})) count[c.range] += c.count;
        double n = 0;
        for (auto& kv : count) n += kv.second;
        row("rainflow counting", "ASTM E1049, Fig. 6", "number of cycles", n, 4.0);
        row("S-N curve, detail category 71", "EN 1993-1-9", "cycles at the fatigue limit", ec3_cycles(std::cbrt(0.4) * 71, 71), 5e6);
    }
    std::fclose(out);
    return 0;
}
