// Finite element verification: (1) tip deflection of a cantilever (400 x 20 x 40 mm, end load)
// against Timoshenko beam theory under mesh refinement for hex8 and hex20; (2) iterations and
// time of the conjugate gradient solver with algebraic multigrid and with Jacobi preconditioning
// against the number of unknowns. Writes fem_convergence.csv and fem_scaling.csv.
// Copyright (c) 2026 Karol Tucki. MIT License.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>

#include "tillbridge/fem.hpp"

using namespace tb;

static Solution cantilever(int order, int n, bool amg, double& tip, double& sec) {
    const Vec3 L{400.0, 20.0, 40.0};
    Mesh m = box_mesh(L, 10 * n, n, 2 * n, {}, order);
    std::vector<int64_t> fixed, end;
    for (size_t i = 0; i < m.n_nodes(); i++) {
        if (m.X[i].x < 1e-9) fixed.push_back(int64_t(i));
        if (m.X[i].x > L.x - 1e-9) end.push_back(int64_t(i));
    }
    std::map<int64_t, Vec3> loads;
    for (auto i : end) loads[i] = {0, 0, -1000.0 / end.size()};
    SolverOptions o;
    o.amg = amg; o.rtol = 1e-10; o.maxiter = 100000;
    auto t0 = std::chrono::steady_clock::now();
    Solution s = solve_static(m, 2.1e5, 0.3, fixed, loads, o);
    sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    tip = 0;
    for (auto i : end) tip -= s.U[i].z / end.size();
    return s;
}

int main() {
    const double E = 2.1e5, b = 20, h = 40, L = 400, P = 1000, I = b * h * h * h / 12, G = E / 2.6, A = b * h;
    const double w = P * L * L * L / (3 * E * I) + P * L / (5.0 / 6.0 * G * A);
    std::FILE* f = std::fopen("fem_convergence.csv", "w");
    std::fprintf(f, "element,n,nodes,tip_mm,timoshenko_mm,rel_error\n");
    for (int order : {1, 2})
        for (int n : {1, 2, 3, 4, 6, 8}) {
            if (order == 2 && n > 6) continue;
            double tip, sec;
            Solution s = cantilever(order, n, true, tip, sec);
            size_t nn = s.U.size();
            std::fprintf(f, "%s,%d,%zu,%.6f,%.6f,%.6e\n", order == 1 ? "hex8" : "hex20", n, nn, tip, w, tip / w - 1);
            std::printf("%s n=%d nodes=%zu tip=%.5f error=%.3e\n", order == 1 ? "hex8" : "hex20", n, nn, tip, tip / w - 1);
        }
    std::fclose(f);
    f = std::fopen("fem_scaling.csv", "w");
    std::fprintf(f, "n,dofs,amg_iterations,amg_s,jacobi_iterations,jacobi_s\n");
    for (int n : {1, 2, 3, 4, 5, 6}) {
        double tip, ta, tj;
        Solution a = cantilever(2, n, true, tip, ta);
        Solution j = cantilever(2, n, false, tip, tj);
        std::fprintf(f, "%d,%zu,%d,%.3f,%d,%.3f\n", n, 3 * a.U.size(), a.iterations, ta, j.iterations, tj);
        std::printf("n=%d dofs=%zu AMG %d it %.2f s | Jacobi %d it %.2f s\n", n, 3 * a.U.size(), a.iterations, ta, j.iterations, tj);
    }
    std::fclose(f);
    return 0;
}
