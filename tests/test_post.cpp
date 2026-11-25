// Rainflow counting, S-N curve, Miner sum, load transfer.
// Copyright (c) 2026 Karol Tucki. MIT License.
#include <cmath>
#include <map>

#include "testing.hpp"
#include "tillbridge/post.hpp"

using namespace tb;

TEST(rainflow_astm_example) {
    // ASTM E1049-85, Fig. 6: -2, 1, -3, 5, -1, 3, -4, 4, -2
    std::map<double, double> count;
    for (const Cycle& c : rainflow({-2, 1, -3, 5, -1, 3, -4, 4, -2})) count[c.range] += c.count;
    std::map<double, double> expected{{3, 0.5}, {4, 1.5}, {6, 0.5}, {8, 1.0}, {9, 0.5}};
    CHECK(count == expected);
}

TEST(rainflow_constant_amplitude) {
    std::vector<double> x(4001);
    for (size_t i = 0; i < x.size(); i++) x[i] = std::sin(20 * kPi * i / 4000.0);
    double n = 0;
    for (const Cycle& c : rainflow(x)) if (std::fabs(c.range - 2) < 1e-3) n += c.count;
    CHECK_NEAR(n, 10, 0.5);
}

TEST(sn_curve_knees) {
    double dsc = 71.0;
    CHECK_REL(ec3_cycles(dsc, dsc), 2e6, 1e-12);
    CHECK_REL(ec3_cycles(0.737 * dsc, dsc), 5e6, 1e-3);
    double kl = std::pow(0.4, 1.0 / 3.0) * std::pow(5e6 / 1e8, 1.0 / 5.0);   // 0.737 and 0.405 are rounded
    CHECK_NEAR(kl, 0.405, 5e-4);
    CHECK_REL(ec3_cycles(kl * dsc * 1.0000001, dsc), 1e8, 1e-4);
    CHECK(std::isinf(ec3_cycles(0.40 * dsc, dsc)));
    CHECK(std::isfinite(ec3_cycles(0.40 * dsc, dsc, false)));
}

TEST(miner_life_scaling) {
    // doubling the stress scales the life by 1/8 on the m = 3 branch
    std::vector<double> F(2001);
    for (size_t i = 0; i < F.size(); i++) F[i] = 1000 * (1 + 0.5 * std::sin(2 * kPi * 5 * (10.0 * i / 2000)));
    FatigueResult a = fatigue_life(F, 15.0, 120.0, 71.0, 1.5), b = fatigue_life(F, 15.0, 240.0, 71.0, 1.5);
    CHECK_REL(a.life_km / b.life_km, 8.0, 1e-6);
    CHECK_NEAR(a.cycles, 50, 1);
}

TEST(distribute_reproduces_wrench) {
    Mesh m;
    Rng rng(0);
    for (int i = 0; i < 40; i++) m.X.push_back({rng.normal(), rng.normal(), rng.normal()});
    std::vector<int64_t> nodes;
    for (int i = 10; i < 40; i++) nodes.push_back(i);
    Vec3 F{100, -50, 20}, M{3, 7, -2}, ref{0.5, 0, 1};
    auto f = distribute(m, nodes, F, M, ref);
    Vec3 Fs, Ms;
    for (auto& [i, fi] : f) { Fs += fi; Ms += cross(m.X[i] - ref, fi); }
    CHECK(norm(Fs - F) < 1e-9);
    CHECK(norm(Ms - M) < 1e-9);
}

TEST(wrench_rotation) {
    double R[3][3] = {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}};
    Vec3 F, M;
    wrench({1, 0, 0}, {0, 0, 1}, R, {}, {}, F, M);
    CHECK(norm(F - Vec3{0, 1, 0}) < 1e-15);
    CHECK(norm(M - cross(Vec3{0, 0, 1}, Vec3{0, 1, 0})) < 1e-15);
}

TEST(centre_of_pressure_weighted) {
    Tool t;
    t.tri = {{Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0}}, {Vec3{2, 0, 0}, Vec3{3, 0, 0}, Vec3{2, 1, 0}}};
    std::vector<TriangleLoad> L(2);
    L[0].Fn = 1; L[1].Fn = 3;
    Vec3 c = centre_of_pressure(t, L);
    Vec3 e = (t.centroid(0) * 1 + t.centroid(1) * 3) * 0.25;
    CHECK(norm(c - e) < 1e-12);
}

TEST(facet_wrench_matches_resultant) {
    // two facets: the wrench about a point equals the sum of the facet forces at their centroids
    Tool t;
    t.tri.push_back({Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0}});
    t.tri.push_back({Vec3{0, 0, 1}, Vec3{1, 0, 1}, Vec3{0, 1, 1}});
    std::vector<TriangleLoad> L(2);
    L[0].F = {10, 0, 0}; L[1].F = {-4, 3, 0};
    double R[3][3] = {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}};   // 90 deg about z
    Vec3 F, M;
    facet_wrench(t, L, 0, 2, R, {0, 0, 0}, 1000.0, {0, 0, 0}, F, M);
    CHECK(norm(F - Vec3{-3, 6, 0}) < 1e-12);
    Vec3 c0{-1000.0 / 3, 1000.0 / 3, 0}, c1{-1000.0 / 3, 1000.0 / 3, 1000};
    Vec3 Me = cross(c0, Vec3{0, 10, 0}) + cross(c1, Vec3{-3, -4, 0});
    CHECK(norm(M - Me) < 1e-9);
}
