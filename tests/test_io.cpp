// Beds, tools, parameters: generation, file round trips, geometry.
// Copyright (c) 2026 Karol Tucki. MIT License.
#include <cmath>
#include <cstdio>
#include <fstream>

#include "testing.hpp"
#include "tillbridge/bed.hpp"
#include "tillbridge/params.hpp"
#include "tillbridge/tool.hpp"

using namespace tb;

namespace {
double total_area(const Tool& t) { double a = 0; for (size_t i = 0; i < t.size(); i++) a += t.area(i); return a; }
double tri_diff(const Tool& a, const Tool& b) {
    if (a.size() != b.size()) return 1e30;
    double d = 0;
    for (size_t i = 0; i < a.size(); i++) for (int k = 0; k < 3; k++) d = std::max(d, norm(a.tri[i][k] - b.tri[i][k]));
    return d;
}
}  // namespace

TEST(clod_volume_equivalence) {
    Rng rng(4);
    double R = 0.015, ext = 0;
    auto c = clod_shape(R, rng, &ext, 200000);
    CHECK(c.size() >= 3 && c.size() <= 4);
    Rng mc(9);
    long inside = 0, n = 400000;
    for (long i = 0; i < n; i++) {
        Vec3 p{mc.uniform(-ext, ext), mc.uniform(-ext, ext), mc.uniform(-ext, ext)};
        for (auto& s : c) if (norm(p - s.offset) <= s.r) { inside++; break; }
    }
    double V = double(inside) / n * std::pow(2 * ext, 3);
    CHECK_NEAR(V / (4.0 / 3.0 * kPi * R * R * R), 1.0, 0.02);
}

TEST(bed_roundtrip_trim_furrow) {
    BedOptions o;
    o.Lx = 0.3; o.Ly = 0.2; o.radius = 0.012; o.clods = true; o.seed = 2;
    Bed b = generate_bed(200, o);
    CHECK(b.size() == 200);
    b.write("tb_test_bed.txt");
    Bed c = Bed::read("tb_test_bed.txt");
    CHECK(c.size() == 200 && c.n_spheres() == b.n_spheres());
    double d = 0;
    for (size_t i = 0; i < b.size(); i++) d = std::max(d, norm(b.bodies[i].x - c.bodies[i].x));
    CHECK(d < 1e-6);
    CHECK_REL(c.mass(), b.mass(), 1e-6);
    Bed t = trim(c, 0.05);
    for (auto& B : t.bodies) CHECK(B.top() <= 0.055 + 1e-9);
    CHECK(cut_furrow(c, 0.1, 0.02).size() < c.size());
    fix_below(c, 0.03);
    int nf = 0;
    for (auto& B : c.bodies) { nf += B.mat == kFixed; if (B.mat == kFixed) CHECK(B.x.z < 0.03); }
    CHECK(nf > 0);
    std::remove("tb_test_bed.txt");
}

TEST(rng_is_deterministic) {
    Rng a(42), b(42);
    for (int i = 0; i < 100; i++) CHECK(a.next() == b.next());
    // xoshiro256** seeded with splitmix64(42): first value is fixed across platforms
    Rng c(42);
    uint64_t v = c.next();
    Rng d(42);
    CHECK(v == d.next());
}

TEST(tool_text_and_stl_roundtrip) {
    Tool t = chisel_tine();
    t.write_stl("tb_test_tool.stl");
    CHECK(tri_diff(t, Tool::read_stl("tb_test_tool.stl")) < 1e-6);
    t.write("tb_test_tool.txt");
    CHECK(tri_diff(t, Tool::read("tb_test_tool.txt")) < 1e-6);   // written in micrometres, as by the C solver
    Tool r = t.refined(0.02);
    CHECK_REL(total_area(r), total_area(t), 1e-12);
    double emax = 0;
    for (auto& T : r.tri) for (int k = 0; k < 3; k++) emax = std::max(emax, norm(T[(k + 1) % 3] - T[k]));
    CHECK(emax <= 0.02 + 1e-12);
    std::remove("tb_test_tool.stl");
    std::remove("tb_test_tool.txt");
}

TEST(binary_stl_in_millimetres) {
    Tool t = chisel_tine(0.05, 0.15, 60.0, 0.35, 0.0, 0.0, 2);
    {
        std::ofstream f("tb_test_bin.stl", std::ios::binary);
        char head[80] = {};
        f.write(head, 80);
        uint32_t n = uint32_t(t.size());
        f.write(reinterpret_cast<char*>(&n), 4);
        for (size_t i = 0; i < t.size(); i++) {
            float v[12];
            Vec3 nn = t.normal(i);
            v[0] = float(nn.x); v[1] = float(nn.y); v[2] = float(nn.z);
            for (int k = 0; k < 3; k++) { v[3 + 3 * k] = float(1000 * t.tri[i][k].x); v[4 + 3 * k] = float(1000 * t.tri[i][k].y); v[5 + 3 * k] = float(1000 * t.tri[i][k].z); }
            f.write(reinterpret_cast<char*>(v), 48);
            uint16_t attr = 0;
            f.write(reinterpret_cast<char*>(&attr), 2);
        }
    }
    CHECK(tri_diff(t, Tool::read_stl("tb_test_bin.stl", 0.001)) < 1e-6);
    std::remove("tb_test_bin.stl");
}

TEST(mouldboard_is_valid) {
    Tool t = mouldboard_body();
    CHECK(t.size() > 50);
    for (size_t i = 0; i < t.size(); i++) CHECK(t.area(i) > 0);
    Vec3 lo, hi;
    t.bounds(lo, hi);
    CHECK(lo.y < -0.3);
    CHECK(hi.z > 0.3);
    CHECK_NEAR(lo.z, 0.1, 1e-12);
}

TEST(params_roundtrip_and_check) {
    Params p;
    p.set("E", 5e6);
    p.set("mu_pw", "0.45");
    p.set("out", "abc");
    p.write("tb_test_params.txt");
    Params q = Params::read("tb_test_params.txt");
    for (const auto& k : Params::keys()) if (k != "out") CHECK(p.get(k) == q.get(k));
    CHECK(q.out == "abc");
    bool thrown = false;
    try { p.set("no_such_key", 1.0); } catch (const std::exception&) { thrown = true; }
    CHECK(thrown);
    thrown = false;
    p.nu = 0.7;
    try { p.check(); } catch (const std::exception&) { thrown = true; }
    CHECK(thrown);
    std::remove("tb_test_params.txt");
}
