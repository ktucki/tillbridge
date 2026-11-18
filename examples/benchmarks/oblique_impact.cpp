// Oblique impact of a sphere on a wall over the whole range of incidence (after Chung and Ooi
// 2011, test 6): rebound angular velocity and tangential velocity against rigid-body theory with
// Coulomb friction. Writes oblique_impact.csv (columns: vt/vn, simulated and theoretical
// omega R / vn and vt'/vn).
// Usage: oblique_impact [OUT.csv]
// Copyright (c) 2026 Karol Tucki. MIT License.
#include <cmath>
#include <cstdio>
#include <string>

#include "tillbridge/dem.hpp"

using namespace tb;

int main(int argc, char** argv) {
    const std::string out = argc > 1 ? argv[1] : "oblique_impact.csv";
    const double R = 0.015, rho = 2600.0, m = rho * 4.0 / 3.0 * kPi * R * R * R, mu = 0.3, vn = 1.0;
    std::FILE* f = std::fopen(out.c_str(), "w");
    if (!f) { std::perror(out.c_str()); return 1; }
    std::fprintf(f, "vt_vn,e_sim,wR_vn_sim,wR_vn_theory,vt_out_vn_sim,vt_out_vn_theory\n");
    for (int k = 0; k <= 40; k++) {
        double vt = 0.1 * k;
        Params p;
        p.E = 2.6e8; p.nu = 0.3; p.mu_pp = 0; p.mu_pw = mu; p.e = 0.5; p.mu_r = 0; p.k_coh = 0; p.g = 0;
        p.dt = 1e-6; p.t_end = 0.02; p.t_tool = 99; p.v_tool = 0; p.out_dt = 1e-3; p.Lx = 2; p.Ly = 0.4; p.Lz = 0.4;
        Bed b;
        Body B;
        B.x = {0.5, 0.2, R + 0.002}; B.m = m; B.I = 0.4 * m * R * R; B.spheres = {Sphere{Vec3{}, R}}; B.v = {vt, 0, -vn};
        b.bodies.push_back(B);
        RunResult r = run_dem(p, b, Tool{}, false);
        const Body& o = r.state.bodies[0];
        double e = o.v.z / vn;
        // rigid sphere, Coulomb friction: sliding throughout if vt > 7/2 mu (1 + e) vn, else ends rolling
        double lim = 3.5 * mu * (1 + e) * vn;
        double wR = vt > lim ? 2.5 * mu * (1 + e) * vn : 5.0 / 7.0 * vt;
        double vto = vt > lim ? vt - mu * (1 + e) * vn : 5.0 / 7.0 * vt;
        std::fprintf(f, "%.3f,%.6f,%.6f,%.6f,%.6f,%.6f\n", vt / vn, e, o.w.y * R / vn, wR / vn, o.v.x / vn, vto / vn);
    }
    std::fclose(f);
    std::printf("written %s\n", out.c_str());
    return 0;
}
