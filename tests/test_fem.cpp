// Finite element verification: patch tests, cantilever against beam theory, AMG against Jacobi,
// ANSYS mesh reading.
// Copyright (c) 2026 Karol Tucki. MIT License.
#include <cmath>
#include <cstdio>
#include <fstream>

#include "testing.hpp"
#include "tillbridge/ansys.hpp"
#include "tillbridge/gmsh.hpp"
#include "tillbridge/fem.hpp"

using namespace tb;

namespace {
struct Cantilever {
    Mesh m;
    Solution s;
    std::vector<int64_t> end;
};
Cantilever cantilever(int order, int nx, int ny, int nz, bool amg = true, double P = 1000.0) {
    Cantilever c;
    Vec3 L{400.0, 20.0, 40.0};
    c.m = box_mesh(L, nx, ny, nz, {}, order);
    std::vector<int64_t> fixed;
    for (size_t i = 0; i < c.m.n_nodes(); i++) {
        if (c.m.X[i].x < 1e-9) fixed.push_back(int64_t(i));
        if (c.m.X[i].x > L.x - 1e-9) c.end.push_back(int64_t(i));
    }
    std::map<int64_t, Vec3> loads;
    for (auto i : c.end) loads[i] = {0, 0, -P / c.end.size()};
    SolverOptions o;
    o.amg = amg; o.rtol = 1e-10;
    c.s = solve_static(c.m, 2.1e5, 0.3, fixed, loads, o);
    return c;
}
double tip(const Cantilever& c) { double w = 0; for (auto i : c.end) w -= c.s.U[i].z; return w / c.end.size(); }
}  // namespace

TEST(cantilever_deflection_hex20) {
    // Timoshenko: P L^3 / (3 E I) + P L / (k G A)
    Cantilever c = cantilever(2, 20, 2, 4);
    double E = 2.1e5, b = 20, h = 40, L = 400, P = 1000, I = b * h * h * h / 12, G = E / 2.6, A = b * h;
    double w = P * L * L * L / (3 * E * I) + P * L / (5.0 / 6.0 * G * A);
    CHECK_REL(tip(c), w, 0.02);
    CHECK_NEAR(c.s.reaction.x, 0, 1e-6 * P);
    CHECK_NEAR(c.s.reaction.z, P, 1e-6 * P);
}

TEST(cantilever_deflection_hex8_stiff) {
    // linear hexahedra lock in bending: deflection below the beam value, within 30 %
    Cantilever c = cantilever(1, 20, 2, 4);
    double E = 2.1e5, b = 20, h = 40, L = 400, P = 1000, I = b * h * h * h / 12, G = E / 2.6, A = b * h;
    double w = P * L * L * L / (3 * E * I) + P * L / (5.0 / 6.0 * G * A);
    CHECK(tip(c) < w && tip(c) > 0.7 * w);
}

TEST(bending_stress_hex20) {
    Cantilever c = cantilever(2, 40, 2, 6);
    double P = 1000, sig = P * 200.0 * 20.0 / (20.0 * 40.0 * 40.0 * 40.0 / 12);
    auto& G = c.m.groups[0];
    double best = 0;
    for (size_t e = 0; e < G.size(); e++) {
        Vec3 xc;
        for (int a = 0; a < 8; a++) xc += c.m.X[G.conn[e * 20 + a]];
        xc = xc * (1.0 / 8);
        if (std::fabs(xc.x - 200.0) < 6 && (xc.z < 5 || xc.z > 35)) best = std::max(best, c.s.vm_elem[0][e]);
    }
    CHECK_REL(best, sig, 0.10);
}

TEST(amg_matches_jacobi) {
    Cantilever a = cantilever(2, 30, 3, 6, true), b = cantilever(2, 30, 3, 6, false);
    double dmax = 0, umax = 0;
    for (size_t i = 0; i < a.s.U.size(); i++) { dmax = std::max(dmax, norm(a.s.U[i] - b.s.U[i])); umax = std::max(umax, norm(b.s.U[i])); }
    CHECK(dmax < 1e-6 * umax);
    CHECK(a.s.iterations < 150);
    CHECK(a.s.iterations < b.s.iterations);
}

TEST(patch_test_distorted_hex20) {
    // uniform tension prescribed on the boundary of a distorted mesh is reproduced exactly inside
    Mesh m = box_mesh({10, 10, 10}, 2, 2, 2, {}, 2);
    std::vector<char> inner(m.n_nodes());
    for (size_t i = 0; i < m.n_nodes(); i++) {
        Vec3& x = m.X[i];
        inner[i] = std::fabs(x.x - 5) < 4.9 && std::fabs(x.y - 5) < 4.9 && std::fabs(x.z - 5) < 4.9;
        if (inner[i]) x += Vec3{0.7, -0.4, 0.3};
    }
    double E = 2.1e5, nu = 0.3, eps = 1e-3;
    CsrMatrix K = assemble_stiffness(m, E, nu);
    size_t n = 3 * m.n_nodes();
    std::vector<double> uex(n);
    for (size_t i = 0; i < m.n_nodes(); i++) { uex[3 * i] = eps * m.X[i].x; uex[3 * i + 1] = -nu * eps * m.X[i].y; uex[3 * i + 2] = -nu * eps * m.X[i].z; }
    // residual of the exact field at the inner dofs must vanish
    std::vector<double> r;
    K.multiply(uex, r);
    double rmax = 0, kmax = 0;
    for (size_t i = 0; i < m.n_nodes(); i++)
        for (int c = 0; c < 3; c++) if (inner[i]) rmax = std::max(rmax, std::fabs(r[3 * i + c]));
    for (double v : K.val) kmax = std::max(kmax, std::fabs(v));
    CHECK(rmax < 1e-9 * kmax * eps * 10);
}

TEST(tet10_volume_and_patch) {
    // unit cube split into six tet10 elements
    Vec3 C[8] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    int T4[6][4] = {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6}, {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}};
    Mesh m;
    for (auto& c : C) m.X.push_back(c);
    std::map<std::pair<int, int>, int64_t> mid;
    auto md = [&](int a, int b) {
        auto k = std::minmax(a, b);
        auto it = mid.find(k);
        if (it != mid.end()) return it->second;
        m.X.push_back((C[a] + C[b]) * 0.5);
        return mid[k] = int64_t(m.X.size() - 1);
    };
    ElementGroup g{ElementKind::Tet10, {}};
    for (auto& t : T4) {
        int a = t[0], b = t[1], c = t[2], d = t[3];
        for (int64_t v : {int64_t(a), int64_t(b), int64_t(c), int64_t(d), md(a, b), md(b, c), md(c, a), md(a, d), md(b, d), md(c, d)}) g.conn.push_back(v);
    }
    m.groups.push_back(g);
    CHECK_NEAR(m.volume(), 1.0, 1e-12);
    // uniaxial tension: fix x = 0 face in x, one node fully, load x = 1 face
    std::vector<int64_t> fixed;
    std::map<int64_t, Vec3> loads;
    double E = 1000, nu = 0.25, sig = 1.0;
    // exact solution check through the stiffness residual (as for hex20)
    CsrMatrix K = assemble_stiffness(m, E, nu);
    std::vector<double> u(3 * m.n_nodes()), r;
    for (size_t i = 0; i < m.n_nodes(); i++) { u[3 * i] = sig / E * m.X[i].x; u[3 * i + 1] = -nu * sig / E * m.X[i].y; u[3 * i + 2] = -nu * sig / E * m.X[i].z; }
    K.multiply(u, r);
    // total force on the x = 1 face equals sigma * area
    double fx = 0;
    for (size_t i = 0; i < m.n_nodes(); i++) if (m.X[i].x > 1 - 1e-12) fx += r[3 * i];
    CHECK_NEAR(fx, sig * 1.0, 1e-12);
}

TEST(ansys_cdb_reader) {
    Mesh m = box_mesh({2, 1, 1}, 2, 1, 1, {}, 1);
    {
        std::ofstream f("tb_test_mesh.cdb");
        char buf[512];
        f << "NBLOCK,6,SOLID," << m.n_nodes() << "," << m.n_nodes() << "\n(3i9,6e21.13e3)\n";
        for (size_t i = 0; i < m.n_nodes(); i++) {
            std::snprintf(buf, sizeof buf, "%9zu%9d%9d%21.13e%21.13e%21.13e\n", i + 1, 0, 0, m.X[i].x, m.X[i].y, m.X[i].z);
            f << buf;
        }
        f << "N,R5.3,LOC,     -1,\nEBLOCK,19,SOLID,2,2\n(19i9)\n";
        for (size_t e = 0; e < 2; e++) {
            int head[11] = {1, 1, 1, 1, 0, 0, 0, 0, 8, 0, int(e + 1)};
            for (int v : head) { std::snprintf(buf, sizeof buf, "%9d", v); f << buf; }
            for (int a = 0; a < 8; a++) { std::snprintf(buf, sizeof buf, "%9lld", (long long)m.groups[0].conn[e * 8 + a] + 1); f << buf; }
            f << "\n";
        }
        f << "       -1\nCMBLOCK,FIXED,NODE,4\n(8i10)\n         1        -3         5         7\nMP,EX,1,210000\nMP,NUXY,1,0.3\n";
    }
    AnsysModel a = AnsysModel::read("tb_test_mesh.cdb");
    Mesh b = a.mesh();
    CHECK(b.n_nodes() == m.n_nodes());
    CHECK(b.n_elements() == 2);
    double d = 0;
    for (size_t i = 0; i < m.n_nodes(); i++) d = std::max(d, norm(b.X[i] - m.X[i]));
    CHECK(d < 1e-12);
    CHECK(a.components["FIXED"] == std::vector<long>({1, 2, 3, 5, 7}));
    CHECK_NEAR(a.materials["EX"], 210000.0, 0);
    CHECK_NEAR(b.volume(), 2.0, 1e-12);
    std::remove("tb_test_mesh.cdb");
}

// One 20-node hexahedron and one 10-node tetrahedron written in the Gmsh node order (version 4.1):
// after reading, every mid-side node of the solver order must lie midway between its corner nodes.
TEST(gmsh_reader_node_order) {
    {
        std::ofstream f("tb_test_mesh.msh");
        f << "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$PhysicalNames\n2\n2 7 \"BASE\"\n3 1 \"BODY\"\n$EndPhysicalNames\n";
        f << "$Entities\n0 0 1 2\n5 0 0 0 1 1 0 1 7 0\n1 0 0 0 1 1 1 1 1 0\n2 2 0 0 3 1 1 0 0\n$EndEntities\n";
        // hexahedron corners (Gmsh = VTK order), then Gmsh edge order 01 03 04 12 15 23 26 37 45 47 56 67
        double c[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
        int ge[12][2] = {{0, 1}, {0, 3}, {0, 4}, {1, 2}, {1, 5}, {2, 3}, {2, 6}, {3, 7}, {4, 5}, {4, 7}, {5, 6}, {6, 7}};
        // tetrahedron, Gmsh edges 01 12 20 30 32 31
        double t[4][3] = {{2, 0, 0}, {3, 0, 0}, {2, 1, 0}, {2, 0, 1}};
        int te[6][2] = {{0, 1}, {1, 2}, {2, 0}, {3, 0}, {3, 2}, {3, 1}};
        f << "$Nodes\n2 30 1 30\n3 1 0 20\n";
        for (int i = 1; i <= 20; i++) f << i << "\n";
        for (int i = 0; i < 8; i++) f << c[i][0] << " " << c[i][1] << " " << c[i][2] << "\n";
        for (auto& e : ge) f << 0.5 * (c[e[0]][0] + c[e[1]][0]) << " " << 0.5 * (c[e[0]][1] + c[e[1]][1]) << " " << 0.5 * (c[e[0]][2] + c[e[1]][2]) << "\n";
        f << "3 2 0 10\n";
        for (int i = 21; i <= 30; i++) f << i << "\n";
        for (int i = 0; i < 4; i++) f << t[i][0] << " " << t[i][1] << " " << t[i][2] << "\n";
        for (auto& e : te) f << 0.5 * (t[e[0]][0] + t[e[1]][0]) << " " << 0.5 * (t[e[0]][1] + t[e[1]][1]) << " " << 0.5 * (t[e[0]][2] + t[e[1]][2]) << "\n";
        f << "$EndNodes\n$Elements\n3 3 1 3\n3 1 17 1\n1";
        for (int i = 1; i <= 20; i++) f << " " << i;
        f << "\n3 2 11 1\n2";
        for (int i = 21; i <= 30; i++) f << " " << i;
        f << "\n2 5 3 1\n3 1 2 3 4\n$EndElements\n";
    }
    AnsysModel a = read_gmsh("tb_test_mesh.msh");
    CHECK(a.solids.size() == 2);
    // solver (VTK) edge order
    int he[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    int tv[6][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};
    double d = 0;
    const auto& H = a.solids[0];
    for (int k = 0; k < 12; k++) d = std::max(d, norm(a.nodes[H[8 + k]] - 0.5 * (a.nodes[H[he[k][0]]] + a.nodes[H[he[k][1]]])));
    const auto& T = a.solids[1];
    for (int k = 0; k < 6; k++) d = std::max(d, norm(a.nodes[T[4 + k]] - 0.5 * (a.nodes[T[tv[k][0]]] + a.nodes[T[tv[k][1]]])));
    CHECK(d < 1e-12);
    CHECK(a.components["BASE"] == std::vector<long>({1, 2, 3, 4}));
    CHECK(a.components["BODY"].size() == 20);
    Mesh m = a.mesh();
    CHECK_NEAR(m.volume(), 1.0 + 1.0 / 6.0, 1e-12);
    std::remove("tb_test_mesh.msh");
}

// ANSYS-type model written as .msh (version 2.2) and read back: same nodes, elements, components and
// the same static solution.
TEST(gmsh_write_read_roundtrip) {
    Mesh box = box_mesh({100, 10, 20}, 10, 2, 2, {}, 2);
    AnsysModel a;
    for (size_t i = 0; i < box.n_nodes(); i++) a.nodes[long(i) + 1] = box.X[i];
    const auto& g = box.groups[0];
    for (size_t e = 0; e < g.size(); e++) {
        std::vector<long> conn;
        for (int k = 0; k < 20; k++) conn.push_back(long(g.conn[e * 20 + k]) + 1);
        a.solids.push_back(conn);
    }
    for (size_t i = 0; i < box.n_nodes(); i++) {
        if (box.X[i].x < 1e-9) a.components["CLAMP"].push_back(long(i) + 1);
        if (box.X[i].x > 100 - 1e-9) a.components["TIP"].push_back(long(i) + 1);
    }
    write_gmsh(a, "tb_test_rt.msh");
    AnsysModel b = read_model("tb_test_rt.msh");
    CHECK(b.nodes.size() == a.nodes.size());
    CHECK(b.solids == a.solids);
    CHECK(b.components["CLAMP"] == a.components["CLAMP"]);
    CHECK(b.components["TIP"] == a.components["TIP"]);
    auto solve = [](const AnsysModel& M) {
        Mesh m = M.mesh();
        std::map<int64_t, Vec3> loads;
        auto tip = M.component("TIP");
        for (auto i : tip) loads[i] = {0, 0, -100.0 / tip.size()};
        return solve_static(m, 2.1e5, 0.3, M.component("CLAMP"), loads);
    };
    Solution s1 = solve(a), s2 = solve(b);
    double d = 0;
    for (size_t i = 0; i < s1.U.size(); i++) d = std::max(d, norm(s1.U[i] - s2.U[i]));
    CHECK(d < 1e-12);
    std::remove("tb_test_rt.msh");
}
