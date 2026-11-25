// Regression against the reference C solver (tbdem v1.3): a chisel tine 6 cm deep in a settled
// 10 cm bed of 622 spheres. With the same compiler settings and one thread the C++ solver
// reproduces the reference force history to the printed digits; across compilers (different
// floating-point contraction) the trajectories of a granular flow diverge slowly, so the mean
// forces are compared with a tolerance.
// Copyright (c) 2026 Karol Tucki. MIT License.
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "testing.hpp"
#include "tillbridge/dem.hpp"

#ifndef TB_TEST_DATA
#define TB_TEST_DATA "tests/data"
#endif

using namespace tb;

TEST(regression_chisel_against_reference) {
    const std::string d = std::string(TB_TEST_DATA) + "/regression/";
    Params p = Params::read(d + "params.txt");
    p.out = "tb_test_regr";
#ifdef _OPENMP
    int nt = omp_get_max_threads();
    omp_set_num_threads(1);
#endif
    RunResult r = run_dem(p, Bed::read(d + "bed.txt"), Tool::read(d + "tine.txt"), false);
#ifdef _OPENMP
    omp_set_num_threads(nt);
#endif
    CHECK(!r.unstable);
    std::ifstream f(d + "reference_force.csv");
    std::string line;
    std::getline(f, line);
    std::vector<ForceSample> ref;
    while (std::getline(f, line)) {
        for (char& c : line) if (c == ',') c = ' ';
        std::istringstream is(line);
        ForceSample s;
        is >> s.t >> s.Fx >> s.Fy >> s.Fz >> s.Ekin >> s.x_tool;
        ref.push_back(s);
    }
    CHECK(ref.size() == r.history.size());
    if (ref.size() != r.history.size()) return;
    // mean force in the averaging window
    double a[3] = {0, 0, 0}, b[3] = {0, 0, 0};
    int n = 0, same = 0;
    for (size_t i = 0; i < ref.size(); i++) {
        const ForceSample &x = r.history[i], &y = ref[i];
        if (std::fabs(x.Fx - y.Fx) < 6e-5 && std::fabs(x.Fz - y.Fz) < 6e-5) same++;
        if (y.x_tool < p.avg_x0 || y.x_tool > p.avg_x1) continue;
        a[0] += x.Fx; a[1] += x.Fy; a[2] += x.Fz; b[0] += y.Fx; b[1] += y.Fy; b[2] += y.Fz; n++;
    }
    CHECK(n > 10);
    CHECK_REL(a[0] / n, b[0] / n, 0.05);
    CHECK_REL(a[2] / n, b[2] / n, 0.05);
    std::printf("  %d of %zu force samples identical to the reference (to the printed digits)\n", same, ref.size());
}

TEST(energy_accounting_closes) {
    // the work of the tool equals draft x distance, and the energy taken up by the contacts plus
    // the change of kinetic and potential energy accounts for it within the integration error
    const std::string d = std::string(TB_TEST_DATA) + "/regression/";
    Params p = Params::read(d + "params.txt");
    p.out = "tb_test_energy";
    RunResult r = run_dem(p, Bed::read(d + "bed.txt"), Tool::read(d + "tine.txt"), false);
    CHECK(!r.unstable);
    const EnergyAccount& e = r.energy;
    CHECK_NEAR(e.distance, p.avg_x1 - p.avg_x0, 2 * p.dt * p.v_tool);
    CHECK_REL(e.tool, r.draft() * e.distance, 2e-3);
    CHECK(e.bulk > 0 && e.tool_contacts > 0);
    CHECK(std::fabs(e.closure()) < 0.03 * e.tool);
}
