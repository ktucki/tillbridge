// Contact-law benchmarks with analytical solutions (after Chung and Ooi 2011, Granular Matter 13,
// 643-656) and solver checks. They test the coupling of rotation and sliding at the contact point,
// which a centre-velocity formulation gets wrong.
// Copyright (c) 2026 Karol Tucki. MIT License.
#include <cmath>

#include "testing.hpp"
#include "tillbridge/dem.hpp"

using namespace tb;

namespace {
const double R = 0.015, RHO = 2600.0;
const double M = RHO * 4.0 / 3.0 * kPi * R * R * R, I = 0.4 * M * R * R;

Bed spheres(const std::vector<Vec3>& x, const std::vector<Vec3>& v = {}, const std::vector<Vec3>& w = {}) {
    Bed b;
    for (size_t i = 0; i < x.size(); i++) {
        Body B;
        B.x = x[i]; B.m = M; B.I = I; B.spheres = {Sphere{Vec3{}, R}};
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

RunResult run(const Params& p, const Bed& b, const Tool& t = {}) { return run_dem(p, b, t, false); }
}  // namespace

TEST(static_hertz_overlap) {
    // a sphere resting on the floor: overlap from k_n delta^1.5 = m g with E* = E/(1-nu^2) (rigid floor)
    Params p = base();
    p.e = 0.1; p.t_end = 0.3; p.dt = 2e-6;
    RunResult r = run(p, spheres({{0.5, 0.2, R}}));
    double Ee = effective_modulus(p.E, p.nu, p.E_wall, p.nu_wall);
    double delta = std::pow(M * p.g / (4.0 / 3.0 * Ee * std::sqrt(R)), 2.0 / 3.0);
    CHECK_REL(R - r.state.bodies[0].x.z, delta, 0.01);
}

TEST(normal_restitution_independent_of_speed) {
    // Tsuji et al. (1992) damping of the Hertz model: the rebound ratio does not depend on the
    // impact speed. The damping ratio is computed from e with the linear-model formula, so the
    // realised value deviates for strongly dissipative contacts (0.10 below e = 0.5, 0.06 above).
    for (double e : {0.3, 0.5, 0.8}) {
        double rr[2];
        int k = 0;
        for (double vn : {1.5, 3.0}) {
            Params p = base();
            p.g = 0; p.e = e; p.t_end = 0.02;
            RunResult r = run(p, spheres({{0.5, 0.2, R + 0.001}}, {{0, 0, -vn}}));
            rr[k++] = r.state.bodies[0].v.z / vn;
        }
        CHECK_NEAR(rr[0], rr[1], 0.01);
        CHECK_NEAR(rr[1], e, e < 0.5 ? 0.10 : 0.06);
    }
}

TEST(oblique_impact_sliding_regime) {
    // dv_t = -mu (1 + e) v_n and d(omega) = 5/2 mu (1 + e) v_n / R; stiff material so that overlap << R
    for (double vt : {2.0, 3.0}) {
        double mu = 0.3, vn = 1.0;
        Params p = base();
        p.E = 2.6e8; p.mu_pw = mu; p.e = 0.5; p.g = 0; p.t_end = 0.02;
        RunResult r = run(p, spheres({{0.5, 0.2, R + 0.002}}, {{vt, 0, -vn}}));
        const Body& b = r.state.bodies[0];
        double e = b.v.z / vn;
        CHECK(e > 0.3);
        CHECK_REL(b.v.x, vt - mu * (1 + e) * vn, 2e-3);
        CHECK_REL(std::fabs(b.w.y), 2.5 * mu * (1 + e) * vn / R, 0.015);
        CHECK(b.w.y > 0);
    }
}

TEST(spinning_sphere_set_down) {
    // ends rolling without slipping: v = 2/7 R omega0
    double w0 = 100.0;
    Params p = base();
    p.mu_pw = 0.5; p.e = 0.3; p.t_end = 0.4; p.dt = 2e-6;
    RunResult r = run(p, spheres({{0.5, 0.2, R - 1e-5}}, {}, {{0, w0, 0}}));
    const Body& b = r.state.bodies[0];
    CHECK_REL(b.v.x, 2.0 / 7.0 * R * w0, 0.01);
    CHECK_REL(b.w.y * R, b.v.x, 0.01);
}

TEST(tool_plate_drags_sphere) {
    // a horizontal tool plate moving under a resting sphere: v = 2/7 v_tool
    double vtool = 0.5, z = 0.2;
    Tool t;
    Vec3 P0{-1.0, 0.0, z}, P1{1.5, 0.0, z}, P2{1.5, 0.4, z}, P3{-1.0, 0.4, z};
    t.tri = {{P0, P1, P2}, {P0, P2, P3}};
    t = t.refined(0.05);
    Params p = base();
    p.mu_pw = 0.5; p.e = 0.3; p.t_tool = 0; p.v_tool = vtool; p.t_end = 0.3; p.dt = 2e-6; p.avg_x0 = -1; p.avg_x1 = 2;
    RunResult r = run(p, spheres({{0.5, 0.2, z + R - 4e-5}}), t);
    CHECK_REL(r.state.bodies[0].v.x, 2.0 / 7.0 * vtool, 0.01);
    // energy bookkeeping: the plate does work on the sphere; the friction power on the tool is dissipative
    double Pt = 0;
    for (auto& L : r.triangles) Pt += L.Pt;
    CHECK(Pt >= -1e-9);
}

TEST(no_rolling_resistance_on_walls) {
    double v0 = 0.3;
    Params p = base();
    p.mu_pw = 0.5; p.e = 0.3; p.mu_r = 0.3; p.t_end = 0.5; p.dt = 2e-6;
    RunResult r = run(p, spheres({{0.3, 0.2, R - 1e-5}}, {{v0, 0, 0}}, {{0, v0 / R, 0}}));
    CHECK_REL(r.state.bodies[0].v.x, v0, 0.01);
}

TEST(binary_collision) {
    double vrel = 1.0;
    Params p = base();
    p.mu_pp = 0.5; p.e = 0.6; p.g = 0; p.t_end = 0.02;
    RunResult r = run(p, spheres({{0.5, 0.2, 0.2}, {0.5 + 2 * R + 0.001, 0.2, 0.2}}, {{vrel / 2, 0, 0}, {-vrel / 2, 0, 0}}));
    const Body &a = r.state.bodies[0], &b = r.state.bodies[1];
    CHECK_NEAR(a.v.x + b.v.x, 0.0, 1e-9);
    double e = (b.v.x - a.v.x) / vrel;
    CHECK(e > 0.55 && e < 0.70);
    CHECK_NEAR(std::fabs(a.w.x) + std::fabs(a.w.y) + std::fabs(a.w.z), 0.0, 1e-9);
}

TEST(oblique_binary_collision_sliding) {
    double mu = 0.2, vn = 1.0, vt = 1.6;
    Params p = base();
    p.E = 2.6e11; p.mu_pp = mu; p.e = 0.6; p.g = 0; p.t_end = 0.0008; p.dt = 1e-8;
    RunResult r = run(p, spheres({{0.5, 0.2, 0.2}, {0.5 + 2 * R + 1e-7, 0.2, 0.2}}, {{vn / 2, vt / 2, 0}, {-vn / 2, -vt / 2, 0}}));
    const Body &a = r.state.bodies[0], &b = r.state.bodies[1];
    double Jn = M * (b.v.x - (-vn / 2));
    CHECK_REL(vt / 2 - a.v.y, mu * Jn / M, 0.01);
    double dw = 2.5 * mu * Jn / (M * R);
    CHECK_REL(std::fabs(a.w.z), dw, 0.01);
    CHECK_REL(std::fabs(b.w.z), dw, 0.01);
    CHECK((a.w.z > 0) == (b.w.z > 0));
}

TEST(energy_growth_flags_unstable) {
    Params p = base();
    p.t_end = 0.2; p.dt = 1e-5; p.out_dt = 1e-3; p.ke_growth = 3.0;
    CHECK(run(p, spheres({{0.5, 0.2, 0.3}})).unstable);
    p.ke_growth = 0; p.t_end = 0.05;
    CHECK(!run(p, spheres({{0.5, 0.2, 0.3}})).unstable);
}

TEST(fixed_bodies_do_not_move) {
    Params p = base();
    p.t_end = 0.05; p.dt = 1e-5;
    Bed b = spheres({{0.5, 0.2, 0.1}});
    b.bodies[0].mat = kFixed;
    RunResult r = run(p, b);
    CHECK_NEAR(r.state.bodies[0].x.z, 0.1, 1e-15);
}

TEST(parallel_bond_holds_and_breaks) {
    // two bonded spheres hanging under gravity: a strong bond holds, a weak one breaks
    for (double strength : {1e9, 1.0}) {
        Params p = base();
        p.g = 9.81; p.t_end = 0.02; p.dt = 1e-6; p.bond_E = 1e7; p.bond_sc = strength; p.bond_tc = strength;
        Bed b = spheres({{0.5, 0.2, 0.30}, {0.5, 0.2, 0.30 - 2 * R}});
        b.bodies[0].mat = kFixed;
        RunResult r = run(p, b);
        CHECK(r.bonds.size() == 1);
        if (strength > 1) { CHECK(r.bonds[0].alive); CHECK_NEAR(r.state.bodies[1].x.z, 0.30 - 2 * R, 1e-4); }
        else { CHECK(!r.bonds[0].alive); CHECK(r.state.bodies[1].x.z < 0.30 - 2 * R - 1e-3); }
    }
}

TEST(closest_point_on_triangle_regions) {
    Vec3 a{0, 0, 0}, b{1, 0, 0}, c{0, 1, 0};
    auto q = closest_on_triangle({0.2, 0.2, 1}, a, b, c);
    CHECK_NEAR(norm(q - Vec3{0.2, 0.2, 0}), 0, 1e-15);
    q = closest_on_triangle({-1, -1, 0}, a, b, c);
    CHECK_NEAR(norm(q - a), 0, 1e-15);
    q = closest_on_triangle({2, -0.5, 0}, a, b, c);
    CHECK_NEAR(norm(q - b), 0, 1e-15);
    q = closest_on_triangle({1, 1, 0}, a, b, c);
    CHECK_NEAR(norm(q - Vec3{0.5, 0.5, 0}), 0, 1e-15);
    q = closest_on_triangle({0.5, -1, 0}, a, b, c);
    CHECK_NEAR(norm(q - Vec3{0.5, 0, 0}), 0, 1e-15);
}

TEST(rayleigh_time_value) {
    // T_R = pi R sqrt(rho/G) / (0.1631 nu + 0.8766)
    double G = 2.6e7 / 2.6;
    CHECK_REL(rayleigh_time(0.015, 2600, 2.6e7, 0.3), kPi * 0.015 * std::sqrt(2600 / G) / (0.1631 * 0.3 + 0.8766), 1e-12);
}
