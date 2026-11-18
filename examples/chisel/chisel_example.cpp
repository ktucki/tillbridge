// Worked example: chisel tine in a soil bin, from particles to fatigue life, through the library API.
//
//  1. generate and settle a bed of spheres, level it to 0.20 m;
//  2. run a reference simulation with known parameters and take its draft and vertical force as
//     "measurements" (a synthetic calibration problem with a known answer);
//  3. recover the cohesion and the soil-steel friction by Bayesian calibration (the particle
//     friction is fixed, as when it is taken from an angle-of-repose test: the draft alone
//     cannot separate particle friction from cohesion);
//  4. transfer the calibrated tool load to a finite element model of the tine shank;
//  5. compute the fatigue life of the shank from the force history.
//
// Usage: chisel_example [OUTDIR]   (default ./out; about 15 min on four cores)
// Steps that have finished are not repeated (bed, reference run, calibration log).
// Figures: python make_figures.py OUTDIR
//
// Copyright (c) 2026 Karol Tucki. MIT License.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <sys/stat.h>

#include "tillbridge/calibrate.hpp"
#include "tillbridge/dem.hpp"
#include "tillbridge/fem.hpp"
#include "tillbridge/post.hpp"

#ifdef _WIN32
#include <direct.h>
static void make_dir(const std::string& d) { _mkdir(d.c_str()); }
#else
static void make_dir(const std::string& d) { mkdir(d.c_str(), 0755); }
#endif
static bool exists(const std::string& f) { return std::ifstream(f).good(); }

using namespace tb;

int main(int argc, char** argv) {
    const std::string out = (argc > 1 ? std::string(argv[1]) : std::string("out")) + "/";
    make_dir(out);
    const double Lx = 1.0, Ly = 0.5, depth = 0.20, radius = 0.012, speed = 1.5;
    const double X0 = 0.30, X1 = 0.75;                 // averaging window (steady state, away from the end wall)
    const double true_coh = 0.9e5, true_mupw = 0.5;    // the "unknown" soil

    Params base;
    base.E = 2.6e6; base.nu = 0.3; base.mu_pp = 0.55; base.e = 0.4; base.mu_r = 0.2; base.dt = 5e-5;
    base.Lx = Lx; base.Ly = Ly; base.Lz = 0.35; base.out_dt = 0.002;

    // ---------------------------------------------------------------- 1. bed
    const std::string bed_file = out + "bed.txt";
    if (!exists(bed_file)) {
        BedOptions o;
        o.Lx = Lx; o.Ly = Ly; o.radius = radius; o.seed = 11;
        Bed column = generate_bed(particles_for_depth(Lx, Ly, depth + 0.01, radius), o);
        Params p = base;
        p.mu_pp = 0.5; p.mu_pw = 0.5; p.k_coh = 0; p.t_end = 1.6; p.t_tool = 99; p.out_dt = 0.05; p.out = out + "settle";
        RunResult s = run_dem(p, column, Tool{});
        Bed b = trim(s.state, depth);
        b.write(bed_file);
        std::printf("bed: %zu particles, bulk density %.0f kg/m3\n", b.size(), b.mass() / (Lx * Ly * depth));
    }
    const Bed bed = Bed::read(bed_file);

    // ---------------------------------------------------------------- tool
    // a wide tine (120 mm, 45 deg): several particles across the face keep the run-to-run scatter
    // of the draft near 1 %, and the vertical force responds to the soil-steel friction
    const double rake = 45.0, width = 0.12;
    auto make_tine = [&](double y) { return chisel_tine(width, 0.12, rake, 0.35, y, 0.0).translated({0, 0, depth - 0.30}).refined(0.03); };
    const Tool tine = make_tine(Ly / 2);
    tine.write(out + "tine.txt");
    tine.write_stl(out + "tine.stl");
    const double t_run = 0.1 + 0.8 / speed;    // the tool travels 0.8 m

    auto simulate_with = [&](const Tool& tool, const std::vector<double>& x, const std::string& label) {
        Params p = base;
        p.k_coh = x[0]; p.mu_pw = x[1]; p.t_end = t_run; p.t_tool = 0.1; p.v_tool = speed;
        p.avg_x0 = X0; p.avg_x1 = X1; p.ke_growth = 3.0; p.out = out + label;
        RunResult r = run_dem(p, bed, tool);
        if (r.unstable) return std::map<std::string, double>{{"draft", 0.0}, {"vertical", 0.0}};   // scored as no load
        return std::map<std::string, double>{{"draft", r.draft()}, {"vertical", -r.F_mean.z}};
    };
    auto simulate = [&](const std::vector<double>& x, const std::string& label) { return simulate_with(tine, x, label); };

    // ---------------------------------------------------------------- 2. reference "measurement"
    // As in a soil bin, the measurement is the mean of four passes (here the tine shifted sideways by
    // -15, -5, 5 and 15 mm); a single DEM run scatters by a few per cent in the draft.
    const double shifts[4] = {-0.015, -0.005, 0.005, 0.015};
    double Dr[4], Vr[4];
    for (int k = 0; k < 4; k++) {
        const std::string lab = "reference_" + std::to_string(k + 1);
        if (!exists(out + lab + "_summary.txt")) simulate_with(make_tine(Ly / 2 + shifts[k]), {true_coh, true_mupw}, lab);
        RunFiles f = RunFiles::read(out + lab);
        Dr[k] = -f.F_mean().x; Vr[k] = -f.F_mean().z;
    }
    double D_ref = 0, V_ref = 0, sD = 0, sV = 0;
    for (int k = 0; k < 4; k++) { D_ref += Dr[k] / 4; V_ref += Vr[k] / 4; }
    for (int k = 0; k < 4; k++) { sD += (Dr[k] - D_ref) * (Dr[k] - D_ref) / 3; sV += (Vr[k] - V_ref) * (Vr[k] - V_ref) / 3; }
    sD = std::sqrt(sD); sV = std::sqrt(sV);
    std::printf("reference (mean of 4 passes): draft %.1f N (sd %.1f), vertical %.1f N (sd %.1f)\n", D_ref, sD, V_ref, sV);

    // ---------------------------------------------------------------- 3. calibration
    Calibration cal(simulate, {"k_coh", "mu_pw"}, {0.0, 0.2}, {1.8e5, 0.9});
    cal.set_targets({{"draft", {D_ref, 1.0}}, {"vertical", {V_ref, 1.0}}});
    cal.set_log(out + "calibration.log");
    cal.set_batch(2);
    cal.set_seed(3);
    cal.run(8, 6);
    // Scattered results: the estimate is the minimum of the expected misfit of the surrogates, not
    // the best single run.
    Evaluation est = cal.estimate();
    std::printf("estimate from draft and vertical force: k_coh %.0f Pa, mu_pw %.3f\n", est.x[0], est.x[1]);
    // Region of parameters consistent with the measurement: surrogate means of the draft and the
    // vertical force against the reference means, with the scatter of the mean of four passes.
    // With the soil-steel friction measured separately (sliding test, mu_pw = 0.5), the cohesion
    // follows from the draft along that line, with an interval from the same scatter.
    const double eD = sD / 2, eV = sV / 2, mu_measured = 0.5;
    double c_cal = 0, c_lo = 1e30, c_hi = -1, r_best = 1e30;
    {
        std::vector<std::vector<double>> G;
        for (int i = 0; i <= 72; i++) for (int j = 0; j <= 70; j++) G.push_back({i * 2500.0, 0.2 + j * 0.01});
        auto P = cal.predict(G);
        std::FILE* fm = std::fopen((out + "calibration_map.csv").c_str(), "w");
        std::fprintf(fm, "k_coh,mu_pw,draft,draft_sd,vertical,vertical_sd,chi2\n");
        for (size_t k = 0; k < G.size(); k++) {
            auto d = P["draft"][k], v = P["vertical"][k];
            double c2 = (d.first - D_ref) * (d.first - D_ref) / (eD * eD + d.second * d.second) +
                        (v.first - V_ref) * (v.first - V_ref) / (eV * eV + v.second * v.second);
            std::fprintf(fm, "%.1f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", G[k][0], G[k][1], d.first, d.second, v.first, v.second, c2);
            if (std::fabs(G[k][1] - mu_measured) < 1e-9) {
                double r = std::fabs(d.first - D_ref);
                if (r < r_best) { r_best = r; c_cal = G[k][0]; }
                if (r < std::sqrt(eD * eD + d.second * d.second)) { c_lo = std::min(c_lo, G[k][0]); c_hi = std::max(c_hi, G[k][0]); }
            }
        }
        std::fclose(fm);
    }
    Evaluation best;
    best.x = {c_cal, mu_measured};
    if (!exists(out + "calibrated_summary.txt")) simulate(best.x, "calibrated");   // run with the calibrated parameters
    best.label = "calibrated";
    std::printf("calibrated with mu_pw = %.2f: k_coh %.0f Pa (%.0f-%.0f, true %.0f)\n", mu_measured, c_cal, c_lo, c_hi, true_coh);


    // ---------------------------------------------------------------- 4. load transfer and FE of the shank
    RunFiles res = RunFiles::read(out + best.label);
    Vec3 cop = centre_of_pressure(tine, res.triangles);
    // FE model of the tine: a plate of 20 mm x 120 mm along the tine face, 350 mm long, clamped at the
    // upper end (frame mounting); FE frame in mm, z along the tine, origin at the tine tip.
    const double L_sh = 350.0, t_sh = 20.0, w_sh = 1000.0 * width;
    Mesh mesh = box_mesh({t_sh, w_sh, L_sh}, 4, 12, 35, {0.0, -w_sh / 2, 0.0}, 2);
    const double th = rake * kPi / 180.0;
    Vec3 ez{-std::cos(th), 0.0, std::sin(th)}, ey{0.0, 1.0, 0.0}, ex = cross(ey, ez);
    double R[3][3] = {{ex.x, ex.y, ex.z}, {ey.x, ey.y, ey.z}, {ez.x, ez.y, ez.z}};   // FE <- DEM rotation
    Vec3 tip = tine.tri[0][0];
    for (auto& T : tine.tri) for (auto& p : T) if (p.z < tip.z) tip = p;
    Vec3 Rtip{dot(ex, tip), dot(ey, tip), dot(ez, tip)};
    Vec3 F, M;
    wrench(res.F_mean(), cop * 1000.0, R, Rtip * -1000.0, {}, F, M);   // N and N mm about the tip
    std::vector<int64_t> top, loaded;
    for (size_t i = 0; i < mesh.n_nodes(); i++) {
        const Vec3& x = mesh.X[i];
        if (x.z > L_sh - 1e-9) top.push_back(int64_t(i));
        if (x.x < 1e-9 && x.z <= 150.0 + 1e-9) loaded.push_back(int64_t(i));   // soil acts on the lower 150 mm of the face
    }
    Solution sol = solve_static(mesh, 2.1e5, 0.3, top, distribute(mesh, loaded, F, M, {}));
    const ElementGroup& G = mesh.groups[0];
    double hot = 0, umax = 0;
    for (size_t e = 0; e < G.size(); e++) {
        double zc = 0;
        for (int a = 0; a < 8; a++) zc += mesh.X[G.conn[e * 20 + a]].z / 8;
        if (zc > L_sh - 20.0) hot = std::max(hot, sol.vm_elem[0][e]);   // structural hot spot at the clamp
    }
    for (auto& u : sol.U) umax = std::max(umax, norm(u));
    write_vtk(out + "shank.vtk", mesh, sol);
    std::printf("FE: %zu nodes, max |u| %.3f mm, hot-spot stress %.1f MPa, AMG-PCG %d iterations\n", mesh.n_nodes(), umax, hot,
                sol.iterations);

    // ---------------------------------------------------------------- 5. fatigue
    FatigueResult life = fatigue_life(force_magnitude(res.history, X0, X1), X1 - X0, hot, 71.0, speed);
    std::printf("fatigue: %.1f cycles/m, max range %.1f MPa, equivalent range %.1f MPa, life %.3g km (%.3g h)\n", life.cycles_per_m,
                life.max_range, life.equivalent_range, life.life_km, life.life_h);

    std::FILE* f = std::fopen((out + "summary.txt").c_str(), "w");
    std::fprintf(f, "draft_reference_sd = %.6g\nvertical_reference_sd = %.6g\n", sD, sV);
    for (int k = 0; k < 4; k++) std::fprintf(f, "draft_pass_%d = %.6g\nvertical_pass_%d = %.6g\n", k + 1, Dr[k], k + 1, Vr[k]);
    std::fprintf(f, "draft_reference = %.6g\nvertical_reference = %.6g\nk_coh_true = %.6g\nmu_pw_true = %.6g\n", D_ref, V_ref,
                 true_coh, true_mupw);
    std::fprintf(f, "k_coh_estimate = %.6g\nmu_pw_estimate = %.6g\n", est.x[0], est.x[1]);
    std::fprintf(f, "k_coh_calibrated = %.6g\nk_coh_low = %.6g\nk_coh_high = %.6g\nmu_pw_measured = %.6g\nbest_run = %s\nruns = %zu\n", c_cal, c_lo,
                 c_hi, mu_measured, best.label.c_str(), cal.data().size());
    std::fprintf(f, "Fx = %.6g\nFy = %.6g\nFz = %.6g\nMx = %.6g\nMy = %.6g\nMz = %.6g\n", F.x, F.y, F.z, M.x, M.y, M.z);
    std::fprintf(f, "cop_x = %.6g\ncop_y = %.6g\ncop_z = %.6g\n", cop.x, cop.y, cop.z);
    std::fprintf(f, "fe_nodes = %zu\npcg_iterations = %d\numax_mm = %.6g\nhot_spot_MPa = %.6g\n", mesh.n_nodes(), sol.iterations,
                 umax, hot);
    std::fprintf(f, "cycles_per_m = %.6g\nmax_range_MPa = %.6g\nequivalent_range_MPa = %.6g\nlife_km = %.6g\nlife_h = %.6g\n",
                 life.cycles_per_m, life.max_range, life.equivalent_range, life.life_km, life.life_h);
    std::fprintf(f, "life_km_no_cutoff = %.6g\n", life.life_km_no_cutoff);
    std::fclose(f);
    return 0;
}
