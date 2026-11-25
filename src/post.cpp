// Copyright (c) 2026 Karol Tucki. MIT License.
#include "tillbridge/post.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>

#include "tillbridge/textio.hpp"

namespace tb {

// ------------------------------------------------------------------ reading results
static bool exists(const std::string& p) { std::ifstream f(p); return bool(f); }

RunFiles RunFiles::read(const std::string& prefix) {
    RunFiles r;
    {
        TextReader t(prefix + "_force.csv");
        std::string s;
        t.line(s);   // header
        while (t.line(s)) {
            auto v = TextReader::parse(s);
            if (v.size() >= 6) r.history.push_back({v[0], v[1], v[2], v[3], v[4], v[5]});
        }
    }
    if (exists(prefix + "_tri.csv")) {
        TextReader t(prefix + "_tri.csv");
        std::string s;
        t.line(s);
        while (t.line(s)) {
            auto v = TextReader::parse(s);
            if (v.size() < 3) continue;
            TriangleLoad L;
            L.Fn = v[1]; L.contacts = v[2];
            if (v.size() >= 6) L.F = {v[3], v[4], v[5]};
            r.triangles.push_back(L);
        }
    }
    if (exists(prefix + "_triW.csv")) {
        TextReader t(prefix + "_triW.csv");
        std::string s;
        t.line(s);
        size_t k = 0;
        while (t.line(s) && k < r.triangles.size()) {
            auto v = TextReader::parse(s);
            if (v.size() >= 3) { r.triangles[k].Pt = v[1]; r.triangles[k].Pn = v[2]; }
            k++;
        }
    }
    std::ifstream f(prefix + "_summary.txt");
    if (!f) throw std::runtime_error("cannot open " + prefix + "_summary.txt");
    std::string line;
    while (std::getline(f, line)) {
        auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string k = line.substr(0, eq);
        k.erase(k.find_last_not_of(" \t") + 1);
        try { r.summary[k] = std::stod(line.substr(eq + 1)); } catch (...) {}
    }
    return r;
}

Vec3 RunFiles::F_mean() const {
    auto g = [&](const char* k) { auto it = summary.find(k); return it == summary.end() ? 0.0 : it->second; };
    return {g("Fx_mean"), g("Fy_mean"), g("Fz_mean")};
}

// ------------------------------------------------------------------ energy and wear
EnergyBalance energy_balance(const RunFiles& r, const Tool& tool, double v) {
    if (v <= 0) throw std::invalid_argument("travel speed must be positive");
    if (r.triangles.size() != tool.size()) throw std::invalid_argument("triangle count of the run and the tool differ");
    EnergyBalance e;
    e.draft = -r.F_mean().x;
    double Pt = 0, Pn = 0, Ftx = 0;
    for (size_t i = 0; i < tool.size(); i++) {
        const auto& L = r.triangles[i];
        Pt += L.Pt; Pn += L.Pn;
        // tangential part of the force on the triangle, projected on the direction of travel
        Vec3 n = tool.normal(i);
        Vec3 Ft = L.F - n * dot(L.F, n);
        Ftx += -Ft.x;
    }
    e.friction = Pt / v;
    e.normal_damping = Pn / v;
    e.soil_and_slice = e.draft - e.friction - e.normal_damping;
    e.tangential_share = e.draft != 0 ? Ftx / e.draft : 0;
    return e;
}

SliceEnergy slice_energy(const std::string& prefix, const Bed& initial, double x0, double x1, double g) {
    if (x1 <= x0) throw std::invalid_argument("empty travel window");
    TextReader t(prefix + "_depart.txt");
    std::string s;
    SliceEnergy e;
    std::vector<char> seen(initial.size(), 0);
    while (t.line(s)) {
        auto v = TextReader::parse(s);
        if (v.size() < 9) continue;
        long b = long(v[1]);
        if (b < 0 || size_t(b) >= initial.size()) throw std::runtime_error("departure record of an unknown body");
        if (v[8] < x0 || v[8] > x1 || seen[b]) continue;   // one record per body (multi-sphere bodies)
        seen[b] = 1;
        const Body& B = initial.bodies[b];
        e.kinetic += 0.5 * B.m * (v[5] * v[5] + v[6] * v[6] + v[7] * v[7]);
        e.lifting += B.m * g * (v[4] - B.x.z);
        e.particles++;
    }
    e.kinetic /= (x1 - x0);
    e.lifting /= (x1 - x0);
    return e;
}

std::vector<double> smooth_on_tool(const Tool& tool, const std::vector<double>& f, int passes) {
    // neighbours: triangles sharing a vertex (vertices matched to 1e-9 m)
    std::map<std::array<long long, 3>, std::vector<size_t>> vtx;
    auto key = [](const Vec3& p) { return std::array<long long, 3>{std::llround(p.x * 1e9), std::llround(p.y * 1e9), std::llround(p.z * 1e9)}; };
    for (size_t i = 0; i < tool.size(); i++) for (auto& p : tool.tri[i]) vtx[key(p)].push_back(i);
    std::vector<std::vector<size_t>> nb(tool.size());
    for (auto& kv : vtx) for (auto i : kv.second) for (auto j : kv.second) nb[i].push_back(j);
    for (auto& v : nb) { std::sort(v.begin(), v.end()); v.erase(std::unique(v.begin(), v.end()), v.end()); }
    std::vector<double> a(tool.size());
    for (size_t i = 0; i < tool.size(); i++) a[i] = tool.area(i);
    std::vector<double> g = f;
    for (int p = 0; p < passes; p++) {
        std::vector<double> h(g.size());
        for (size_t i = 0; i < g.size(); i++) {
            double s = 0, w = 0;
            for (auto j : nb[i]) { s += g[j] * a[j]; w += a[j]; }
            h[i] = w > 0 ? s / w : g[i];
        }
        g.swap(h);
    }
    return g;
}

WearMap wear_map(const RunFiles& r, const Tool& tool, double p_min) {
    if (r.triangles.size() != tool.size()) throw std::invalid_argument("triangle count of the run and the tool differ");
    WearMap w;
    size_t n = tool.size();
    w.pressure.assign(n, 0); w.q_f.assign(n, 0); w.slip_speed.assign(n, 0);
    double pa = 0, aa = 0;
    std::vector<double> vs;
    for (size_t i = 0; i < n; i++) {
        double A = std::max(tool.area(i), 1e-12);
        w.pressure[i] = std::max(r.triangles[i].Fn, 0.0) / A;
        w.q_f[i] = std::max(r.triangles[i].Pt, 0.0) / A;
        w.total_frictional_power += std::max(r.triangles[i].Pt, 0.0);
        // slip speed: frictional power over the tangential force (exact for any mix of sticking and sliding)
        Vec3 n = tool.normal(i);
        Vec3 Ft = r.triangles[i].F - n * dot(r.triangles[i].F, n);
        double ftm = norm(Ft);
        if (w.pressure[i] > p_min && ftm > 0) {
            w.slip_speed[i] = std::max(r.triangles[i].Pt, 0.0) / ftm;
            vs.push_back(w.slip_speed[i]);
            pa += w.pressure[i] * A; aa += A;
        }
    }
    w.mean_pressure = aa > 0 ? pa / aa : 0;
    if (!vs.empty()) { std::nth_element(vs.begin(), vs.begin() + vs.size() / 2, vs.end()); w.median_slip_speed = vs[vs.size() / 2]; }
    return w;
}

// ------------------------------------------------------------------ load transfer
Vec3 centre_of_pressure(const Tool& tool, const std::vector<TriangleLoad>& L) {
    if (L.size() != tool.size()) throw std::invalid_argument("triangle count of the run and the tool differ");
    Vec3 c;
    double w = 0;
    for (size_t i = 0; i < tool.size(); i++) { double f = std::max(L[i].Fn, 0.0); c += tool.centroid(i) * f; w += f; }
    if (w <= 0) throw std::runtime_error("no contact force on the tool");
    return c * (1.0 / w);
}

void wrench(const Vec3& Fd, const Vec3& cop, const double R[3][3], const Vec3& t, const Vec3& ref, Vec3& F, Vec3& M) {
    auto rot = [&](const Vec3& p) {
        return Vec3{R[0][0] * p.x + R[0][1] * p.y + R[0][2] * p.z, R[1][0] * p.x + R[1][1] * p.y + R[1][2] * p.z,
                    R[2][0] * p.x + R[2][1] * p.y + R[2][2] * p.z};
    };
    F = rot(Fd);
    Vec3 p = rot(cop) + t;
    M = cross(p - ref, F);
}

void facet_wrench(const Tool& tool, const std::vector<TriangleLoad>& L, size_t first, size_t end, const double R[3][3],
                  const Vec3& t, double scale, const Vec3& ref, Vec3& F, Vec3& M) {
    if (L.size() != tool.size()) throw std::invalid_argument("triangle count of the run and the tool differ");
    end = std::min(end, tool.size());
    auto rot = [&](const Vec3& p) {
        return Vec3{R[0][0] * p.x + R[0][1] * p.y + R[0][2] * p.z, R[1][0] * p.x + R[1][1] * p.y + R[1][2] * p.z,
                    R[2][0] * p.x + R[2][1] * p.y + R[2][2] * p.z};
    };
    F = Vec3{}; M = Vec3{};
    for (size_t k = first; k < end; k++) {
        Vec3 c = (rot(tool.centroid(k)) + t) * scale, Fk = rot(L[k].F);
        F += Fk; M += cross(c - ref, Fk);
    }
}

std::map<int64_t, Vec3> distribute(const Mesh& mesh, const std::vector<int64_t>& nodes, const Vec3& F, const Vec3& M0, const Vec3& ref) {
    if (nodes.size() < 3) throw std::invalid_argument("load distribution needs at least three nodes");
    Vec3 c;
    for (auto i : nodes) c += mesh.X.at(i);
    c = c * (1.0 / nodes.size());
    Vec3 M = M0 + cross(ref - c, F);
    // S = sum(r.r) I - sum r r^T ; w = S^-1 M ; f_i = F/n + w x r_i
    double S[3][3] = {}, rr = 0;
    for (auto i : nodes) {
        Vec3 r = mesh.X[i] - c;
        rr += dot(r, r);
        for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) S[a][b] -= r[a] * r[b];
    }
    for (int a = 0; a < 3; a++) S[a][a] += rr;
    double det = S[0][0] * (S[1][1] * S[2][2] - S[1][2] * S[2][1]) - S[0][1] * (S[1][0] * S[2][2] - S[1][2] * S[2][0]) +
                 S[0][2] * (S[1][0] * S[2][1] - S[1][1] * S[2][0]);
    if (std::fabs(det) < 1e-30) throw std::runtime_error("node set is collinear: the moment cannot be distributed");
    double Si[3][3] = {{(S[1][1] * S[2][2] - S[1][2] * S[2][1]) / det, (S[0][2] * S[2][1] - S[0][1] * S[2][2]) / det, (S[0][1] * S[1][2] - S[0][2] * S[1][1]) / det},
                       {(S[1][2] * S[2][0] - S[1][0] * S[2][2]) / det, (S[0][0] * S[2][2] - S[0][2] * S[2][0]) / det, (S[0][2] * S[1][0] - S[0][0] * S[1][2]) / det},
                       {(S[1][0] * S[2][1] - S[1][1] * S[2][0]) / det, (S[0][1] * S[2][0] - S[0][0] * S[2][1]) / det, (S[0][0] * S[1][1] - S[0][1] * S[1][0]) / det}};
    Vec3 w{Si[0][0] * M.x + Si[0][1] * M.y + Si[0][2] * M.z, Si[1][0] * M.x + Si[1][1] * M.y + Si[1][2] * M.z,
           Si[2][0] * M.x + Si[2][1] * M.y + Si[2][2] * M.z};
    std::map<int64_t, Vec3> out;
    for (auto i : nodes) out[i] += F * (1.0 / nodes.size()) + cross(w, mesh.X[i] - c);
    return out;
}

std::vector<double> force_magnitude(const std::vector<ForceSample>& h, double x0, double x1) {
    std::vector<double> f;
    for (auto& s : h) if (s.x_tool >= x0 && s.x_tool <= x1) f.push_back(std::sqrt(s.Fx * s.Fx + s.Fy * s.Fy + s.Fz * s.Fz));
    return f;
}

// ------------------------------------------------------------------ fatigue
std::vector<double> turning_points(const std::vector<double>& x) {
    if (x.size() < 3) return x;
    std::vector<double> out{x[0]};
    int last = 0;
    for (size_t i = 0; i + 1 < x.size(); i++) {
        double d = x[i + 1] - x[i];
        int s = d > 0 ? 1 : (d < 0 ? -1 : 0);
        if (s == 0) continue;
        if (last != 0 && s != last) out.push_back(x[i]);
        last = s;
    }
    out.push_back(x.back());
    return out;
}

std::vector<Cycle> rainflow(const std::vector<double>& x) {
    std::vector<double> s;
    std::vector<Cycle> out;
    for (double p : turning_points(x)) {
        s.push_back(p);
        while (s.size() >= 3) {
            size_t n = s.size();
            double X = std::fabs(s[n - 1] - s[n - 2]), Y = std::fabs(s[n - 2] - s[n - 3]);
            if (X < Y) break;
            if (n == 3) {
                out.push_back({Y, (s[1] + s[0]) / 2, 0.5});
                s.erase(s.begin());
            } else {
                out.push_back({Y, (s[n - 2] + s[n - 3]) / 2, 1.0});
                double last = s.back();
                s.resize(n - 3);
                s.push_back(last);
            }
        }
    }
    for (size_t i = 0; i + 1 < s.size(); i++) out.push_back({std::fabs(s[i + 1] - s[i]), (s[i + 1] + s[i]) / 2, 0.5});
    return out;
}

static const double KD = std::cbrt(2.0 / 5.0);                 // dsD / dsc = 0.737
static const double KL = KD * std::pow(5.0 / 100.0, 0.2);       // dsL / dsc = 0.405

double ec3_cycles(double ds, double dsc, bool cutoff) {
    const double inf = std::numeric_limits<double>::infinity();
    double dsd = KD * dsc, dsl = KL * dsc;
    if (ds <= 0) return inf;
    if (cutoff && ds <= dsl) return inf;
    if (ds >= dsd) return 2e6 * std::pow(dsc / ds, 3);
    return 5e6 * std::pow(dsd / ds, 5);
}

double miner(const std::vector<Cycle>& c, double dsc, bool cutoff) {
    double D = 0;
    for (auto& y : c) D += y.count / ec3_cycles(y.range, dsc, cutoff);
    return D;
}

double equivalent_range(const std::vector<Cycle>& c, double m) {
    double s = 0, n = 0;
    for (auto& y : c) { s += y.count * std::pow(y.range, m); n += y.count; }
    return n > 0 ? std::pow(s / n, 1 / m) : 0.0;
}

FatigueResult fatigue_life(const std::vector<double>& F, double distance, double sigma_ref, double dsc, double speed, double F_ref) {
    if (F.empty()) throw std::invalid_argument("empty force history");
    if (distance <= 0) throw std::invalid_argument("distance must be positive");
    if (F_ref <= 0) { for (double f : F) F_ref += f; F_ref /= F.size(); }
    std::vector<double> sig(F.size());
    for (size_t i = 0; i < F.size(); i++) sig[i] = sigma_ref * F[i] / F_ref;
    auto cyc = rainflow(sig);
    FatigueResult r;
    const double inf = std::numeric_limits<double>::infinity();
    for (auto& c : cyc) { r.cycles += c.count; r.max_range = std::max(r.max_range, c.range); }
    r.cycles_per_m = r.cycles / distance;
    r.equivalent_range = equivalent_range(cyc);
    r.ratio_to_cutoff = r.max_range / (KL * dsc);
    r.damage_per_m = miner(cyc, dsc, true) / distance;
    r.damage_per_m_no_cutoff = miner(cyc, dsc, false) / distance;
    r.life_km = r.damage_per_m > 0 ? 1 / r.damage_per_m / 1000 : inf;
    r.life_h = r.damage_per_m > 0 ? 1 / (r.damage_per_m * speed * 3600) : inf;
    r.life_km_no_cutoff = r.damage_per_m_no_cutoff > 0 ? 1 / r.damage_per_m_no_cutoff / 1000 : inf;
    return r;
}

}  // namespace tb
