// tillbridge - command line interface.
// Copyright (c) 2026 Karol Tucki. MIT License.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "tillbridge/ansys.hpp"
#include "tillbridge/gmsh.hpp"
#include "tillbridge/bed.hpp"
#include "tillbridge/calibrate.hpp"
#include "tillbridge/dem.hpp"
#include "tillbridge/fem.hpp"
#include "tillbridge/post.hpp"
#include "tillbridge/tool.hpp"

using namespace tb;

static const char* kUsage = R"(tillbridge - calibrated DEM loads on tillage tools, carried to FE stresses and fatigue life

usage: tillbridge <command> [arguments]

  bed OUT --Lx X --Ly Y --depth D [--radius R] [--clods] [--density RHO] [--bulk RHO_B] [--seed N]
        particles in columns above the box floor, to be settled with `run`
  prepare STATE OUT --depth D [--furrow-edge Y --furrow-bottom Z [--furrow-angle DEG]] [--fix-below Z]
        trim a settled bed, cut the previous furrow, anchor the subsoil
  tool STL OUT [--unit U] [--translate X Y Z] [--max-edge L]
        STL (ASCII or binary) to the solver's triangle file
  shape mouldboard|chisel OUT [--stl FILE]
        parametric demonstration tools
  run BED TOOL|none OUT [--params FILE] [--set KEY=VALUE ...] [--threads N] [--verbose]
        DEM run (TOOL none: settling); outputs OUT_force.csv, _tri.csv, _triW.csv, _state.txt, _vel.txt, _depart.txt, _summary.txt
  energy PREFIX TOOL --speed V [--bed INITIAL --x0 X0 --x1 X1]
        energy balance of the draft work (friction, normal damping, soil; kinetic and lifting work)
  wear PREFIX TOOL [--smooth N] [--out CSV] [--p-min PA]
        pressure, slip speed and frictional power density (wear index) on the tool
  loads PREFIX TOOL [--ref X Y Z] [--facets FIRST END] [--rotation R11..R33 --translation T --scale S]
        mean force and moment of the facet forces about a point and the centre of pressure, optionally
        for a group of facets and in the frame of an FE model (p' = S (R p + T), F' = R F)
  fem MODEL --fixed COMP --loaded COMP --force FX FY FZ [--moment MX MY MZ] [--ref X Y Z]
        [--E E] [--nu NU] [--vtk FILE] [--ccx FILE.inp] [--jacobi]
        static FE analysis of a mesh from ANSYS (.cdb, Mechanical .dat) or Gmsh (.msh), load as a
        distributed coupling; --ccx also writes the same model as a CalculiX input file
  convert MODEL OUT.msh
        mesh and node components of an ANSYS model to the open Gmsh format
  fatigue PREFIX --x0 X0 --x1 X1 --sigma-ref S --speed V [--detail DSC]
        rainflow and Palmgren-Miner fatigue life (EN 1993-1-9 S-N curve)
  calibrate CONFIG
        Bayesian calibration of DEM parameters (see docs/calibration.md)
  version
)";

// ------------------------------------------------------------------ argument parsing
struct Args {
    std::vector<std::string> pos;
    std::map<std::string, std::vector<std::string>> opt;
    Args(int argc, char** argv, int first, const std::map<std::string, int>& nvals) {
        for (int i = first; i < argc; i++) {
            std::string a = argv[i];
            if (a.rfind("--", 0) == 0) {
                auto it = nvals.find(a);
                if (it == nvals.end()) throw std::invalid_argument("unknown option " + a);
                std::vector<std::string> v;
                if (it->second < 0) { while (i + 1 < argc && std::string(argv[i + 1]).rfind("--", 0) != 0) v.push_back(argv[++i]); }
                else for (int k = 0; k < it->second; k++) { if (++i >= argc) throw std::invalid_argument("missing value for " + a); v.push_back(argv[i]); }
                if (it->second < 0) { auto& o = opt[a]; o.insert(o.end(), v.begin(), v.end()); }   // repeatable
                else opt[a] = v;
            } else pos.push_back(a);
        }
    }
    bool has(const std::string& k) const { return opt.count(k) > 0; }
    double num(const std::string& k, int i = 0) const {
        auto it = opt.find(k);
        if (it == opt.end()) throw std::invalid_argument("option " + k + " is required");
        return std::stod(it->second.at(i));
    }
    double num(const std::string& k, double def, int i) const { return has(k) ? num(k, i) : def; }
    std::string str(const std::string& k, const std::string& def = "") const { auto it = opt.find(k); return it == opt.end() ? def : it->second.at(0); }
    Vec3 vec(const std::string& k, const Vec3& def = {}) const { return has(k) ? Vec3{num(k, 0), num(k, 1), num(k, 2)} : def; }
    void need(size_t n, const char* what) const { if (pos.size() != n) throw std::invalid_argument(std::string("usage: tillbridge ") + what); }
};

static std::string fmt(double v) {
    if (std::isinf(v)) return "\"inf\"";
    char b[64];
    std::snprintf(b, sizeof b, "%.6g", v);
    return b;
}
static void print(const std::vector<std::pair<std::string, double>>& kv) {
    std::printf("{\n");
    for (size_t i = 0; i < kv.size(); i++) std::printf("  \"%s\": %s%s\n", kv[i].first.c_str(), fmt(kv[i].second).c_str(), i + 1 < kv.size() ? "," : "");
    std::printf("}\n");
}

// ------------------------------------------------------------------ commands
static int cmd_bed(int argc, char** argv) {
    Args a(argc, argv, 2, {{"--Lx", 1}, {"--Ly", 1}, {"--depth", 1}, {"--radius", 1}, {"--clods", 0}, {"--density", 1}, {"--bulk", 1}, {"--seed", 1}});
    a.need(1, "bed OUT --Lx X --Ly Y --depth D");
    BedOptions o;
    o.Lx = a.num("--Lx"); o.Ly = a.num("--Ly"); o.radius = a.num("--radius", 0.015, 0); o.clods = a.has("--clods");
    o.density = a.num("--density", 2600.0, 0); o.seed = uint64_t(a.num("--seed", 1, 0));
    int n = particles_for_depth(o.Lx, o.Ly, a.num("--depth"), o.radius, o.density, a.num("--bulk", 1500.0, 0));
    Bed b = generate_bed(n, o);
    b.write(a.pos[0]);
    double top = 0;
    for (auto& x : b.bodies) top = std::max(top, x.top());
    std::printf("%zu bodies, %zu spheres, top of the columns %.3f m\n", b.size(), b.n_spheres(), top);
    return 0;
}

static int cmd_prepare(int argc, char** argv) {
    Args a(argc, argv, 2, {{"--depth", 1}, {"--furrow-edge", 1}, {"--furrow-bottom", 1}, {"--furrow-angle", 1}, {"--fix-below", 1}});
    a.need(2, "prepare STATE OUT --depth D");
    Bed b = trim(Bed::read(a.pos[0]), a.num("--depth"));
    if (a.has("--furrow-edge")) b = cut_furrow(b, a.num("--furrow-edge"), a.num("--furrow-bottom"), a.num("--furrow-angle", 45.0, 0));
    if (a.has("--fix-below")) fix_below(b, a.num("--fix-below"));
    b.write(a.pos[1]);
    std::printf("%zu bodies kept\n", b.size());
    return 0;
}

static int cmd_tool(int argc, char** argv) {
    Args a(argc, argv, 2, {{"--unit", 1}, {"--translate", 3}, {"--max-edge", 1}});
    a.need(2, "tool STL OUT");
    Tool t = Tool::read_stl(a.pos[0], a.num("--unit", 1.0, 0)).translated(a.vec("--translate"));
    if (a.has("--max-edge")) t = t.refined(a.num("--max-edge"));
    t.write(a.pos[1]);
    Vec3 lo, hi;
    t.bounds(lo, hi);
    std::printf("%zu triangles, bounds (%.3f %.3f %.3f) .. (%.3f %.3f %.3f) m\n", t.size(), lo.x, lo.y, lo.z, hi.x, hi.y, hi.z);
    return 0;
}

static int cmd_shape(int argc, char** argv) {
    Args a(argc, argv, 2, {{"--stl", 1}});
    a.need(2, "shape mouldboard|chisel OUT");
    Tool t;
    if (a.pos[0] == "mouldboard") t = mouldboard_body();
    else if (a.pos[0] == "chisel") t = chisel_tine();
    else throw std::invalid_argument("unknown shape " + a.pos[0]);
    t.write(a.pos[1]);
    if (a.has("--stl")) t.write_stl(a.str("--stl"), a.pos[0]);
    std::printf("%zu triangles\n", t.size());
    return 0;
}

static int cmd_run(int argc, char** argv) {
    Args a(argc, argv, 2, {{"--params", 1}, {"--set", -1}, {"--threads", 1}, {"--verbose", 0}});
    a.need(3, "run BED TOOL OUT [--params FILE] [--set KEY=VALUE ...]");
    Params p = a.has("--params") ? Params::read(a.str("--params")) : Params{};
    if (a.has("--set"))
        for (auto& kv : a.opt.at("--set")) {
            auto eq = kv.find('=');
            if (eq == std::string::npos) throw std::invalid_argument("--set needs KEY=VALUE, got " + kv);
            p.set(kv.substr(0, eq), kv.substr(eq + 1));
        }
    p.out = a.pos[2];
#ifdef _OPENMP
    if (a.has("--threads")) omp_set_num_threads(int(a.num("--threads")));
#endif
    p.write(p.out + "_params.txt");
    Tool tool = a.pos[1] == "none" ? Tool{} : Tool::read(a.pos[1]);   // "none": settling without a tool
    RunResult r = run_dem(p, Bed::read(a.pos[0]), tool, true, a.has("--verbose"));
    std::vector<std::pair<std::string, double>> kv;
    for (auto& s : r.summary()) kv.push_back(s);
    print(kv);
    return r.unstable ? 3 : 0;
}

static int cmd_energy(int argc, char** argv) {
    Args a(argc, argv, 2, {{"--speed", 1}, {"--bed", 1}, {"--x0", 1}, {"--x1", 1}, {"--g", 1}});
    a.need(2, "energy PREFIX TOOL --speed V");
    RunFiles r = RunFiles::read(a.pos[0]);
    Tool t = Tool::read(a.pos[1]);
    EnergyBalance e = energy_balance(r, t, a.num("--speed"));
    std::vector<std::pair<std::string, double>> kv = {{"draft_J_per_m", e.draft}, {"friction_J_per_m", e.friction},
                                                      {"normal_damping_J_per_m", e.normal_damping}, {"soil_and_slice_J_per_m", e.soil_and_slice},
                                                      {"tangential_share_of_draft", e.tangential_share}};
    if (a.has("--bed")) {
        SliceEnergy s = slice_energy(a.pos[0], Bed::read(a.str("--bed")), a.num("--x0"), a.num("--x1"), a.num("--g", 9.81, 0));
        kv.push_back({"kinetic_J_per_m", s.kinetic});
        kv.push_back({"lifting_J_per_m", s.lifting});
        kv.push_back({"soil_dissipation_J_per_m", e.soil_and_slice - s.kinetic - s.lifting});
        kv.push_back({"departed_particles", double(s.particles)});
    }
    // direct energy accounting of the run (version 2.1 and later)
    auto S = [&](const char* k) { auto it = r.summary.find(k); return it == r.summary.end() ? std::nan("") : it->second; };
    const double L = S("E_distance");
    if (L > 0) {
        kv.push_back({"tool_work_J_per_m", S("E_tool") / L});
        kv.push_back({"tool_contacts_J_per_m", S("E_tool_contacts") / L});
        kv.push_back({"soil_direct_J_per_m", (S("E_bulk") + S("dE_kinetic") + S("dE_potential")) / L});
        kv.push_back({"kinetic_change_J_per_m", S("dE_kinetic") / L});
        kv.push_back({"potential_change_J_per_m", S("dE_potential") / L});
        kv.push_back({"closure_J_per_m", S("E_closure") / L});
        kv.push_back({"closure_relative", S("E_closure") / S("E_tool")});
    }
    print(kv);
    return 0;
}

static int cmd_wear(int argc, char** argv) {
    Args a(argc, argv, 2, {{"--smooth", 1}, {"--out", 1}, {"--p-min", 1}});
    a.need(2, "wear PREFIX TOOL");
    RunFiles r = RunFiles::read(a.pos[0]);
    Tool t = Tool::read(a.pos[1]);
    WearMap w = wear_map(r, t, a.num("--p-min", 1000.0, 0));
    int ns = int(a.num("--smooth", 0, 0));
    if (ns > 0) { w.pressure = smooth_on_tool(t, w.pressure, ns); w.q_f = smooth_on_tool(t, w.q_f, ns); w.slip_speed = smooth_on_tool(t, w.slip_speed, ns); }
    if (a.has("--out")) {
        std::FILE* f = std::fopen(a.str("--out").c_str(), "w");
        if (!f) throw std::runtime_error("cannot write " + a.str("--out"));
        std::fprintf(f, "tri,cx,cy,cz,area,p_Pa,qf_W_m2,vs_m_s\n");
        for (size_t i = 0; i < t.size(); i++) {
            Vec3 c = t.centroid(i);
            std::fprintf(f, "%zu,%.6f,%.6f,%.6f,%.6e,%.6e,%.6e,%.6e\n", i, c.x, c.y, c.z, t.area(i), w.pressure[i], w.q_f[i], w.slip_speed[i]);
        }
        std::fclose(f);
    }
    double qmax = 0;
    for (double q : w.q_f) qmax = std::max(qmax, q);
    print({{"total_frictional_power_W", w.total_frictional_power}, {"mean_pressure_Pa", w.mean_pressure},
           {"median_slip_speed_m_per_s", w.median_slip_speed}, {"max_frictional_power_density_W_per_m2", qmax}});
    return 0;
}

static int cmd_loads(int argc, char** argv) {
    Args a(argc, argv, 2, {{"--ref", 3}, {"--facets", 2}, {"--rotation", 9}, {"--translation", 3}, {"--scale", 1}});
    a.need(2, "loads PREFIX TOOL [--facets FIRST END] [--rotation R11..R33 --translation TX TY TZ --scale S] [--ref X Y Z]");
    RunFiles r = RunFiles::read(a.pos[0]);
    Tool t = Tool::read(a.pos[1]);
    if (r.triangles.size() != t.size()) throw std::invalid_argument("the run and the tool have different numbers of facets");
    size_t f0 = 0, f1 = t.size();
    if (a.has("--facets")) { f0 = size_t(a.num("--facets", 0)); f1 = std::min(t.size(), size_t(a.num("--facets", 1))); }
    // target frame: p' = s (R p + t), F' = R F (e.g. DEM metres to an FE model in millimetres)
    double R[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    if (a.has("--rotation")) for (int i = 0; i < 9; i++) R[i / 3][i % 3] = a.num("--rotation", i);
    const Vec3 tr = a.vec("--translation");
    const double sc = a.num("--scale", 1.0, 0);
    auto P = [&](const Vec3& p) { return Vec3{R[0][0] * p.x + R[0][1] * p.y + R[0][2] * p.z + tr.x, R[1][0] * p.x + R[1][1] * p.y + R[1][2] * p.z + tr.y,
                                              R[2][0] * p.x + R[2][1] * p.y + R[2][2] * p.z + tr.z} * sc; };
    auto V = [&](const Vec3& v) { return Vec3{R[0][0] * v.x + R[0][1] * v.y + R[0][2] * v.z, R[1][0] * v.x + R[1][1] * v.y + R[1][2] * v.z,
                                              R[2][0] * v.x + R[2][1] * v.y + R[2][2] * v.z}; };
    const Vec3 ref = a.vec("--ref");
    // exact wrench of the facet forces; centre of pressure weighted by the normal forces
    Vec3 F, M, cop;
    double fn = 0;
    for (size_t k = f0; k < f1; k++) {
        Vec3 c = P(t.centroid(k)), Fk = V(r.triangles[k].F);
        F += Fk; M += cross(c - ref, Fk);
        cop += c * r.triangles[k].Fn; fn += r.triangles[k].Fn;
    }
    if (fn > 0) cop = cop * (1.0 / fn);
    print({{"facets", double(f1 - f0)}, {"Fx", F.x}, {"Fy", F.y}, {"Fz", F.z}, {"Mx", M.x}, {"My", M.y}, {"Mz", M.z},
           {"cop_x", cop.x}, {"cop_y", cop.y}, {"cop_z", cop.z}});
    return 0;
}

static int cmd_fem(int argc, char** argv) {
    Args a(argc, argv, 2, {{"--fixed", 1}, {"--loaded", 1}, {"--force", 3}, {"--moment", 3}, {"--ref", 3}, {"--E", 1}, {"--nu", 1},
                           {"--vtk", 1}, {"--ccx", 1}, {"--jacobi", 0}, {"--verbose", 0}});
    a.need(1, "fem MODEL --fixed COMP --loaded COMP --force FX FY FZ");
    AnsysModel m = read_model(a.pos[0]);
    Mesh mesh = m.mesh();
    double E = a.has("--E") ? a.num("--E") : (m.materials.count("EX") ? m.materials.at("EX") : throw std::invalid_argument("give --E"));
    double nu = a.has("--nu") ? a.num("--nu") : (m.materials.count("NUXY") ? m.materials.at("NUXY") : 0.3);
    auto fixed = m.component(a.str("--fixed"));
    auto loaded = m.component(a.str("--loaded"));
    Vec3 ref = a.vec("--ref");
    if (!a.has("--ref")) { for (auto i : loaded) ref += mesh.X[i]; ref = ref * (1.0 / loaded.size()); }
    auto loads = distribute(mesh, loaded, a.vec("--force"), a.vec("--moment"), ref);
    if (a.has("--ccx")) write_calculix(a.str("--ccx"), mesh, E, nu, fixed, loads);
    SolverOptions o;
    o.amg = !a.has("--jacobi"); o.verbose = a.has("--verbose");
    Solution s = solve_static(mesh, E, nu, fixed, loads, o);
    if (a.has("--vtk")) write_vtk(a.str("--vtk"), mesh, s);
    double vmax = 0, umax = 0;
    for (double v : s.vm_node) vmax = std::max(vmax, v);
    for (auto& u : s.U) umax = std::max(umax, norm(u));
    print({{"nodes", double(mesh.n_nodes())}, {"elements", double(mesh.n_elements())}, {"max_von_mises_nodal", vmax},
           {"max_displacement", umax}, {"reaction_x", s.reaction.x}, {"reaction_y", s.reaction.y}, {"reaction_z", s.reaction.z},
           {"reaction_moment_x", (s.reaction_moment - cross(ref, s.reaction)).x}, {"reaction_moment_y", (s.reaction_moment - cross(ref, s.reaction)).y},
           {"reaction_moment_z", (s.reaction_moment - cross(ref, s.reaction)).z},
           {"pcg_iterations", double(s.iterations)}});
    return 0;
}

static int cmd_convert(int argc, char** argv) {
    Args a(argc, argv, 2, {});
    a.need(2, "convert MODEL OUT.msh");
    AnsysModel m = read_model(a.pos[0]);
    write_gmsh(m, a.pos[1]);
    print({{"nodes", double(m.nodes.size())}, {"elements", double(m.solids.size())}, {"components", double(m.components.size())}});
    return 0;
}

static int cmd_fatigue(int argc, char** argv) {
    Args a(argc, argv, 2, {{"--x0", 1}, {"--x1", 1}, {"--sigma-ref", 1}, {"--detail", 1}, {"--speed", 1}});
    a.need(1, "fatigue PREFIX --x0 X0 --x1 X1 --sigma-ref S --speed V");
    RunFiles r = RunFiles::read(a.pos[0]);
    double x0 = a.num("--x0"), x1 = a.num("--x1");
    FatigueResult f = fatigue_life(force_magnitude(r.history, x0, x1), x1 - x0, a.num("--sigma-ref"), a.num("--detail", 71.0, 0), a.num("--speed"));
    print({{"cycles", f.cycles}, {"cycles_per_m", f.cycles_per_m}, {"max_range_MPa", f.max_range}, {"equivalent_range_MPa", f.equivalent_range},
           {"ratio_to_cutoff", f.ratio_to_cutoff}, {"damage_per_m", f.damage_per_m}, {"life_km", f.life_km}, {"life_h", f.life_h},
           {"damage_per_m_no_cutoff", f.damage_per_m_no_cutoff}, {"life_km_no_cutoff", f.life_km_no_cutoff}});
    return 0;
}

// Configuration (one item per line, '#' comments):
//   bed FILE | tool FILE | params FILE | set KEY VALUE | param NAME LOW HIGH
//   target QUANTITY VALUE [WEIGHT] | objective QUANTITY | constraint QUANTITY LOW HIGH
//   init N | iter N | batch N | seed N | log FILE | out PREFIX
// Quantities: the summary keys of a run (Fx_mean, Fz_mean, ...) and draft = -Fx_mean.
static int cmd_calibrate(int argc, char** argv) {
    if (argc != 3) throw std::invalid_argument("usage: tillbridge calibrate CONFIG");
    std::ifstream f(argv[2]);
    if (!f) throw std::runtime_error(std::string("cannot open ") + argv[2]);
    std::string bedf, toolf, out = "cal", logf = "calibration.log", objective;
    Params base;
    std::vector<std::string> names;
    std::vector<double> lo, hi;
    std::map<std::string, std::pair<double, double>> targets, constraints;
    int n_init = 8, n_iter = 10, batch = 1;
    uint64_t seed = 0;
    std::string line;
    while (std::getline(f, line)) {
        auto c = line.find('#');
        if (c != std::string::npos) line.resize(c);
        std::istringstream is(line);
        std::string k;
        if (!(is >> k)) continue;
        if (k == "bed") is >> bedf;
        else if (k == "tool") is >> toolf;
        else if (k == "params") { std::string p; is >> p; base = Params::read(p); }
        else if (k == "set") { std::string a, b; is >> a >> b; base.set(a, b); }
        else if (k == "param") { std::string n; double a, b; is >> n >> a >> b; base.get(n); names.push_back(n); lo.push_back(a); hi.push_back(b); }
        else if (k == "target") {   // quantity, measured value, optional weight
            std::string q, ws; double v = 0, w = 1; is >> q >> v;
            if (is.fail()) throw std::invalid_argument("bad line in the calibration file: " + line);
            if (is >> ws) w = std::stod(ws);
            is.clear();
            targets[q] = {v, w};
        }
        else if (k == "objective") is >> objective;
        else if (k == "constraint") { std::string q; double a, b; is >> q >> a >> b; constraints[q] = {a, b}; }
        else if (k == "init") is >> n_init;
        else if (k == "iter") is >> n_iter;
        else if (k == "batch") is >> batch;
        else if (k == "seed") is >> seed;
        else if (k == "log") is >> logf;
        else if (k == "out") is >> out;
        else throw std::invalid_argument("unknown item in the calibration file: " + k);
        if (is.fail()) throw std::invalid_argument("bad line in the calibration file: " + line);
    }
    if (bedf.empty() || toolf.empty() || names.empty()) throw std::invalid_argument("the calibration file needs bed, tool and at least one param");
    Bed bed = Bed::read(bedf);
    Tool tool = Tool::read(toolf);
    auto evaluate = [&](const std::vector<double>& x, const std::string& label) {
        Params p = base;
        for (size_t k = 0; k < names.size(); k++) p.set(names[k], x[k]);
        p.out = out + "_" + label;
        RunResult r = run_dem(p, bed, tool, true);
        auto s = r.summary();
        s["draft"] = r.draft();
        if (r.unstable) {   // numerically unstable parameter set: scored as no load, the search moves away from it
            std::fprintf(stderr, "tillbridge: run %s unstable, scored as zero load\n", p.out.c_str());
            for (auto& kv : s) if (kv.first != "unstable") kv.second = 0;
        }
        return s;
    };
    Calibration cal(evaluate, names, lo, hi);
    if (!targets.empty()) cal.set_targets(targets);
    else cal.set_objective(objective, constraints);
    cal.set_log(logf);
    cal.set_batch(batch);
    cal.set_seed(seed);
    Evaluation b = cal.run(n_init, n_iter), e = cal.estimate();
    std::vector<std::pair<std::string, double>> kv;
    // estimate: minimum of the surrogate means (recommended, DEM results scatter); best: best single run
    for (size_t k = 0; k < names.size(); k++) kv.push_back({names[k], e.x[k]});
    kv.push_back({"score_predicted", e.score});
    for (auto& r : e.result) kv.push_back({"predicted_" + r.first, r.second});
    for (size_t k = 0; k < names.size(); k++) kv.push_back({"best_run_" + names[k], b.x[k]});
    kv.push_back({"best_run_score", b.score});
    kv.push_back({"best_run_feasible", b.feasible ? 1.0 : 0.0});
    print(kv);
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) { std::fputs(kUsage, stderr); return 2; }
    std::string c = argv[1];
    try {
        if (c == "bed") return cmd_bed(argc, argv);
        if (c == "prepare") return cmd_prepare(argc, argv);
        if (c == "tool") return cmd_tool(argc, argv);
        if (c == "shape") return cmd_shape(argc, argv);
        if (c == "run") return cmd_run(argc, argv);
        if (c == "energy") return cmd_energy(argc, argv);
        if (c == "wear") return cmd_wear(argc, argv);
        if (c == "loads") return cmd_loads(argc, argv);
        if (c == "fem") return cmd_fem(argc, argv);
        if (c == "fatigue") return cmd_fatigue(argc, argv);
        if (c == "convert") return cmd_convert(argc, argv);
        if (c == "calibrate") return cmd_calibrate(argc, argv);
        if (c == "version" || c == "--version") { std::printf("tillbridge %s\n", kVersion); return 0; }
        if (c == "help" || c == "--help" || c == "-h") { std::fputs(kUsage, stdout); return 0; }
        // compatibility with the earlier solver interface: tillbridge PARAMS BED TOOL
        if (argc == 4) {
            Params p = Params::read(argv[1]);
            RunResult r = run_dem(p, Bed::read(argv[2]), Tool::read(argv[3]));
            return r.unstable ? 3 : 0;
        }
        std::fprintf(stderr, "tillbridge: unknown command %s\n\n%s", c.c_str(), kUsage);
        return 2;
    } catch (std::exception& e) {
        std::fprintf(stderr, "tillbridge: %s\n", e.what());
        return 1;
    }
}
