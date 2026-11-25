// Copyright (c) 2026 Karol Tucki. MIT License.
#include "tillbridge/fem.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <unordered_map>

namespace tb {

int nodes_per_element(ElementKind k) {
    switch (k) {
        case ElementKind::Hex8: return 8;
        case ElementKind::Hex20: case ElementKind::Hex20Reduced: return 20;
        case ElementKind::Tet10: return 10;
    }
    return 0;
}
const char* element_name(ElementKind k) {
    switch (k) {
        case ElementKind::Hex8: return "hex8";
        case ElementKind::Hex20: return "hex20";
        case ElementKind::Hex20Reduced: return "hex20r";
        case ElementKind::Tet10: return "tet10";
    }
    return "?";
}
size_t Mesh::n_elements() const { size_t n = 0; for (auto& g : groups) n += g.size(); return n; }

// ------------------------------------------------------------------ shape functions
static const double kHex20C[20][3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1}, {-1, -1, 1}, {1, -1, 1}, {1, 1, 1},
                                      {-1, 1, 1},   {0, -1, -1}, {1, 0, -1}, {0, 1, -1},  {-1, 0, -1}, {0, -1, 1}, {1, 0, 1},
                                      {0, 1, 1},    {-1, 0, 1},  {-1, -1, 0}, {1, -1, 0}, {1, 1, 0},   {-1, 1, 0}};

void shape_functions(ElementKind k, double xi, double et, double ze, double* N, double (*d)[3]) {
    if (k == ElementKind::Hex8) {
        for (int i = 0; i < 8; i++) {
            double a = kHex20C[i][0], b = kHex20C[i][1], g = kHex20C[i][2];
            N[i] = (1 + a * xi) * (1 + b * et) * (1 + g * ze) / 8;
            d[i][0] = a * (1 + b * et) * (1 + g * ze) / 8;
            d[i][1] = b * (1 + a * xi) * (1 + g * ze) / 8;
            d[i][2] = g * (1 + a * xi) * (1 + b * et) / 8;
        }
        return;
    }
    if (k == ElementKind::Hex20 || k == ElementKind::Hex20Reduced) {
        for (int i = 0; i < 20; i++) {
            double a = kHex20C[i][0], b = kHex20C[i][1], g = kHex20C[i][2];
            if (i < 8) {
                N[i] = (1 + a * xi) * (1 + b * et) * (1 + g * ze) * (a * xi + b * et + g * ze - 2) / 8;
                d[i][0] = a * (1 + b * et) * (1 + g * ze) * (2 * a * xi + b * et + g * ze - 1) / 8;
                d[i][1] = b * (1 + a * xi) * (1 + g * ze) * (a * xi + 2 * b * et + g * ze - 1) / 8;
                d[i][2] = g * (1 + a * xi) * (1 + b * et) * (a * xi + b * et + 2 * g * ze - 1) / 8;
            } else if (a == 0) {
                N[i] = (1 - xi * xi) * (1 + b * et) * (1 + g * ze) / 4;
                d[i][0] = -2 * xi * (1 + b * et) * (1 + g * ze) / 4;
                d[i][1] = b * (1 - xi * xi) * (1 + g * ze) / 4;
                d[i][2] = g * (1 - xi * xi) * (1 + b * et) / 4;
            } else if (b == 0) {
                N[i] = (1 - et * et) * (1 + a * xi) * (1 + g * ze) / 4;
                d[i][0] = a * (1 - et * et) * (1 + g * ze) / 4;
                d[i][1] = -2 * et * (1 + a * xi) * (1 + g * ze) / 4;
                d[i][2] = g * (1 - et * et) * (1 + a * xi) / 4;
            } else {
                N[i] = (1 - ze * ze) * (1 + a * xi) * (1 + b * et) / 4;
                d[i][0] = a * (1 - ze * ze) * (1 + b * et) / 4;
                d[i][1] = b * (1 - ze * ze) * (1 + a * xi) / 4;
                d[i][2] = -2 * ze * (1 + a * xi) * (1 + b * et) / 4;
            }
        }
        return;
    }
    // Tet10: volume coordinates L0 = 1 - r - s - t, L1 = r, L2 = s, L3 = t
    double L[4] = {1 - xi - et - ze, xi, et, ze};
    const double dL[4][3] = {{-1, -1, -1}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (int i = 0; i < 4; i++) {
        N[i] = L[i] * (2 * L[i] - 1);
        for (int c = 0; c < 3; c++) d[i][c] = (4 * L[i] - 1) * dL[i][c];
    }
    const int e[6][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};
    for (int k2 = 0; k2 < 6; k2++) {
        int a = e[k2][0], b = e[k2][1];
        N[4 + k2] = 4 * L[a] * L[b];
        for (int c = 0; c < 3; c++) d[4 + k2][c] = 4 * (L[a] * dL[b][c] + L[b] * dL[a][c]);
    }
}

struct GaussRule {
    std::vector<std::array<double, 3>> p;
    std::vector<double> w;
};
static GaussRule rule(ElementKind k) {
    GaussRule r;
    if (k == ElementKind::Tet10) {
        double a = 0.5854101966249685, b = 0.1381966011250105;
        r.p = {{b, b, b}, {a, b, b}, {b, a, b}, {b, b, a}};
        r.w.assign(4, 1.0 / 24);
        return r;
    }
    std::vector<double> g, w;
    if (k == ElementKind::Hex20) { g = {-std::sqrt(0.6), 0.0, std::sqrt(0.6)}; w = {5.0 / 9, 8.0 / 9, 5.0 / 9}; }
    else { g = {-1 / std::sqrt(3.0), 1 / std::sqrt(3.0)}; w = {1.0, 1.0}; }
    for (size_t i = 0; i < g.size(); i++)
        for (size_t j = 0; j < g.size(); j++)
            for (size_t l = 0; l < g.size(); l++) { r.p.push_back({g[i], g[j], g[l]}); r.w.push_back(w[i] * w[j] * w[l]); }
    return r;
}

static void elasticity(double E, double nu, double D[6][6]) {
    double l = E * nu / ((1 + nu) * (1 - 2 * nu)), m = E / (2 * (1 + nu));
    for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) D[i][j] = 0;
    for (int i = 0; i < 3; i++) { for (int j = 0; j < 3; j++) D[i][j] = l; D[i][i] += 2 * m; D[3 + i][3 + i] = m; }
}

// Global derivatives at a Gauss point; returns detJ (throws for a non-positive Jacobian).
static double gradients(const Mesh& M, const int64_t* conn, int nn, const double (*dn)[3], double (*dN)[3]) {
    double J[3][3] = {};
    for (int a = 0; a < nn; a++) {
        const Vec3& x = M.X[conn[a]];
        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) J[i][j] += x[i] * dn[a][j];
    }
    double det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
                 J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
    if (det <= 0) throw std::runtime_error("element with a non-positive Jacobian (node order?)");
    double Ji[3][3] = {{(J[1][1] * J[2][2] - J[1][2] * J[2][1]) / det, (J[0][2] * J[2][1] - J[0][1] * J[2][2]) / det, (J[0][1] * J[1][2] - J[0][2] * J[1][1]) / det},
                       {(J[1][2] * J[2][0] - J[1][0] * J[2][2]) / det, (J[0][0] * J[2][2] - J[0][2] * J[2][0]) / det, (J[0][2] * J[1][0] - J[0][0] * J[1][2]) / det},
                       {(J[1][0] * J[2][1] - J[1][1] * J[2][0]) / det, (J[0][1] * J[2][0] - J[0][0] * J[2][1]) / det, (J[0][0] * J[1][1] - J[0][1] * J[1][0]) / det}};
    // dN/dx_i = sum_j dn/dxi_j * (J^-1)_{j i}, with J_{ij} = dx_i/dxi_j
    for (int a = 0; a < nn; a++)
        for (int i = 0; i < 3; i++) dN[a][i] = dn[a][0] * Ji[0][i] + dn[a][1] * Ji[1][i] + dn[a][2] * Ji[2][i];
    return det;
}

// B^T D B contribution: K[3a+i][3b+j]
static void element_stiffness(const Mesh& M, const int64_t* conn, ElementKind k, const double D[6][6], const GaussRule& g,
                              std::vector<double>& Ke) {
    int nn = nodes_per_element(k), nd = 3 * nn;
    Ke.assign(size_t(nd) * nd, 0.0);
    double N[20], dn[20][3], dN[20][3];
    std::vector<double> B(6 * size_t(nd)), DB(6 * size_t(nd));
    for (size_t q = 0; q < g.p.size(); q++) {
        shape_functions(k, g.p[q][0], g.p[q][1], g.p[q][2], N, dn);
        double wd = gradients(M, conn, nn, dn, dN) * g.w[q];
        std::fill(B.begin(), B.end(), 0.0);
        for (int a = 0; a < nn; a++) {
            double x = dN[a][0], y = dN[a][1], z = dN[a][2];
            B[0 * nd + 3 * a] = x; B[1 * nd + 3 * a + 1] = y; B[2 * nd + 3 * a + 2] = z;
            B[3 * nd + 3 * a] = y; B[3 * nd + 3 * a + 1] = x;
            B[4 * nd + 3 * a + 1] = z; B[4 * nd + 3 * a + 2] = y;
            B[5 * nd + 3 * a] = z; B[5 * nd + 3 * a + 2] = x;
        }
        for (int r = 0; r < 6; r++)
            for (int c = 0; c < nd; c++) {
                double s = 0;
                for (int l = 0; l < 6; l++) s += D[r][l] * B[l * nd + c];
                DB[r * nd + c] = s;
            }
        for (int i = 0; i < nd; i++)
            for (int j = 0; j < nd; j++) {
                double s = 0;
                for (int r = 0; r < 6; r++) s += B[r * nd + i] * DB[r * nd + j];
                Ke[size_t(i) * nd + j] += s * wd;
            }
    }
}

double Mesh::volume() const {
    double v = 0, N[20], dn[20][3], dN[20][3];
    for (auto& G : groups) {
        GaussRule g = rule(G.kind);
        int nn = nodes_per_element(G.kind);
        for (size_t e = 0; e < G.size(); e++)
            for (size_t q = 0; q < g.p.size(); q++) {
                shape_functions(G.kind, g.p[q][0], g.p[q][1], g.p[q][2], N, dn);
                v += gradients(*this, &G.conn[e * nn], nn, dn, dN) * g.w[q];
            }
    }
    return v;
}

CsrMatrix assemble_stiffness(const Mesh& M, double E, double nu) {
    const int64_t nn = int64_t(M.n_nodes());
    // node graph
    std::vector<std::vector<int64_t>> adj(nn);
    for (auto& G : M.groups) {
        int k = nodes_per_element(G.kind);
        for (size_t e = 0; e < G.size(); e++)
            for (int a = 0; a < k; a++)
                for (int b = 0; b < k; b++) adj[G.conn[e * k + a]].push_back(G.conn[e * k + b]);
    }
    for (auto& v : adj) { std::sort(v.begin(), v.end()); v.erase(std::unique(v.begin(), v.end()), v.end()); }
    CsrMatrix K;
    K.rows = K.cols = 3 * nn;
    K.ptr.assign(K.rows + 1, 0);
    for (int64_t i = 0; i < nn; i++)
        for (int c = 0; c < 3; c++) K.ptr[3 * i + c + 1] = K.ptr[3 * i + c] + 3 * int64_t(adj[i].size());
    K.idx.resize(K.ptr.back()); K.val.assign(K.ptr.back(), 0.0);
    for (int64_t i = 0; i < nn; i++)
        for (int c = 0; c < 3; c++) {
            int64_t p = K.ptr[3 * i + c];
            for (auto j : adj[i]) for (int d = 0; d < 3; d++) K.idx[p++] = 3 * j + d;
        }
    double D[6][6];
    elasticity(E, nu, D);
    for (auto& G : M.groups) {
        GaussRule g = rule(G.kind);
        int k = nodes_per_element(G.kind), nd = 3 * k;
        const int64_t ne = int64_t(G.size());
#pragma omp parallel
        {
            std::vector<double> Ke;
#pragma omp for schedule(dynamic, 64)
            for (int64_t e = 0; e < ne; e++) {
                const int64_t* conn = &G.conn[e * k];
                element_stiffness(M, conn, G.kind, D, g, Ke);
                for (int a = 0; a < k; a++) {
                    const auto& row = adj[conn[a]];
                    for (int b = 0; b < k; b++) {
                        int64_t pos = int64_t(std::lower_bound(row.begin(), row.end(), conn[b]) - row.begin());
                        for (int i = 0; i < 3; i++) {
                            int64_t base = K.ptr[3 * conn[a] + i] + 3 * pos;
                            for (int j = 0; j < 3; j++) {
#pragma omp atomic
                                K.val[base + j] += Ke[size_t(3 * a + i) * nd + 3 * b + j];
                            }
                        }
                    }
                }
            }
        }
    }
    return K;
}

// ------------------------------------------------------------------ smoothed-aggregation AMG
namespace {
struct SaAmg {
    const CsrMatrix* K = nullptr;
    std::vector<double> Dinv;
    CsrMatrix P, PT;
    DenseCholesky coarse;
    double wj = 0;
    int nsm = 2;

    void setup(const CsrMatrix& Kf, const std::vector<Vec3>& Xf, double h, bool verbose) {
        K = &Kf;
        const int64_t n = int64_t(Xf.size());
        Vec3 lo{1e300, 1e300, 1e300};
        for (auto& x : Xf) for (int c = 0; c < 3; c++) lo[c] = std::min(lo[c], x[c]);
        std::vector<int64_t> agg(n);
        int64_t na = 0;
        for (int attempt = 0; attempt < 20; attempt++) {
            std::unordered_map<int64_t, int64_t> id;
            for (int64_t i = 0; i < n; i++) {
                int64_t a = int64_t(std::floor((Xf[i].x - lo.x) / h)), b = int64_t(std::floor((Xf[i].y - lo.y) / h)),
                        c = int64_t(std::floor((Xf[i].z - lo.z) / h));
                int64_t key = (a * 1000003 + b) * 1000003 + c;
                auto it = id.find(key);
                if (it == id.end()) it = id.emplace(key, int64_t(id.size())).first;
                agg[i] = it->second;
            }
            na = int64_t(id.size());
            if (6 * na <= 6000) break;   // dense coarse problem stays small
            h *= 1.26;
        }
        std::vector<Vec3> cen(na);
        std::vector<int> cnt(na, 0);
        for (int64_t i = 0; i < n; i++) { cen[agg[i]] += Xf[i]; cnt[agg[i]]++; }
        for (int64_t a = 0; a < na; a++) cen[a] = cen[a] * (1.0 / cnt[a]);
        std::vector<Triplet> t;
        for (int64_t i = 0; i < n; i++) {
            Vec3 r = (Xf[i] - cen[agg[i]]) * (1.0 / h);
            for (int d = 0; d < 3; d++) t.push_back({3 * i + d, 6 * agg[i] + d, 1.0});
            const int pairs[3][2] = {{1, 2}, {2, 0}, {0, 1}};
            for (int k = 0; k < 3; k++) {
                int a = pairs[k][0], b = pairs[k][1];
                t.push_back({3 * i + a, 6 * agg[i] + 3 + k, -r[b]});
                t.push_back({3 * i + b, 6 * agg[i] + 3 + k, r[a]});
            }
        }
        CsrMatrix P0 = CsrMatrix::from_triplets(3 * n, 6 * na, t);
        auto dg = Kf.diagonal();
        Dinv.resize(dg.size());
        for (size_t i = 0; i < dg.size(); i++) Dinv[i] = 1.0 / dg[i];
        // spectral radius of D^-1 K by power iteration
        std::vector<double> x(3 * n), y;
        uint64_t s = 12345;
        for (auto& v : x) { s = s * 6364136223846793005ULL + 1442695040888963407ULL; v = double(s >> 11) * 0x1.0p-53; }
        double lam = 1;
        for (int it = 0; it < 15; it++) {
            Kf.multiply(x, y);
            for (size_t i = 0; i < y.size(); i++) y[i] *= Dinv[i];
            lam = norm2(y) / norm2(x);
            double ny = norm2(y);
            for (size_t i = 0; i < y.size(); i++) x[i] = y[i] / ny;
        }
        double omega = 4.0 / 3.0 / lam;
        // P = (I - omega D^-1 K) P0
        CsrMatrix KP = multiply(Kf, P0);
        for (int64_t i = 0; i < KP.rows; i++)
            for (int64_t k = KP.ptr[i]; k < KP.ptr[i + 1]; k++) KP.val[k] *= -omega * Dinv[i];
        std::vector<Triplet> tp;
        for (int64_t i = 0; i < P0.rows; i++) {
            for (int64_t k = P0.ptr[i]; k < P0.ptr[i + 1]; k++) tp.push_back({i, P0.idx[k], P0.val[k]});
            for (int64_t k = KP.ptr[i]; k < KP.ptr[i + 1]; k++) tp.push_back({i, KP.idx[k], KP.val[k]});
        }
        P = CsrMatrix::from_triplets(P0.rows, P0.cols, tp);
        PT = P.transpose();
        CsrMatrix Ac = multiply(PT, multiply(Kf, P));
        int64_t m = Ac.rows;
        std::vector<double> A(size_t(m) * m, 0.0);
        double dmax = 0;
        for (int64_t i = 0; i < m; i++)
            for (int64_t k = Ac.ptr[i]; k < Ac.ptr[i + 1]; k++) { A[i * m + Ac.idx[k]] = Ac.val[k]; if (Ac.idx[k] == i) dmax = std::max(dmax, Ac.val[k]); }
        for (int64_t i = 0; i < m; i++) {
            double d = A[i * m + i];
            A[i * m + i] += 1e-9 * d + (d <= 1e-10 * dmax ? dmax : 0.0);
        }
        coarse.factor(std::move(A), m);
        wj = 0.75 * omega;
        if (verbose) std::fprintf(stderr, "SA-AMG: %lld nodes, %lld aggregates, coarse size %lld\n", (long long)n, (long long)na, (long long)m);
    }

    void apply(const std::vector<double>& b, std::vector<double>& x) const {
        const size_t n = b.size();
        std::vector<double> r(n), Kx, rc;
        x.resize(n);
        for (size_t i = 0; i < n; i++) x[i] = wj * Dinv[i] * b[i];
        auto smooth = [&]() {
            K->multiply(x, Kx);
            for (size_t i = 0; i < n; i++) x[i] += wj * Dinv[i] * (b[i] - Kx[i]);
        };
        for (int k = 0; k < nsm - 1; k++) smooth();
        K->multiply(x, Kx);
        for (size_t i = 0; i < n; i++) r[i] = b[i] - Kx[i];
        PT.multiply(r, rc);
        coarse.solve(rc);
        std::vector<double> e;
        P.multiply(rc, e);
        for (size_t i = 0; i < n; i++) x[i] += e[i];
        for (int k = 0; k < nsm; k++) smooth();
    }
};
}  // namespace

Solution solve_static(const Mesh& M, double E, double nu, const std::vector<int64_t>& fixed, const std::map<int64_t, Vec3>& loads,
                      const SolverOptions& o) {
    auto t0 = std::chrono::steady_clock::now();
    const int64_t nn = int64_t(M.n_nodes());
    CsrMatrix K = assemble_stiffness(M, E, nu);
    std::vector<double> F(3 * nn, 0.0);
    for (auto& l : loads) {
        if (l.first < 0 || l.first >= nn) throw std::out_of_range("load on a node outside the mesh");
        for (int c = 0; c < 3; c++) F[3 * l.first + c] += l.second[c];
    }
    std::vector<char> isfix(nn, 0);
    for (auto i : fixed) { if (i < 0 || i >= nn) throw std::out_of_range("fixed node outside the mesh"); isfix[i] = 1; }
    std::vector<int64_t> free, freenodes;
    for (int64_t i = 0; i < nn; i++)
        if (!isfix[i]) { freenodes.push_back(i); for (int c = 0; c < 3; c++) free.push_back(3 * i + c); }
    if (freenodes.size() == size_t(nn)) throw std::invalid_argument("no fixed nodes: the model can move as a rigid body");
    CsrMatrix Kf = K.principal(free);
    std::vector<double> Ff(free.size()), uf(free.size(), 0.0);
    for (size_t i = 0; i < free.size(); i++) Ff[i] = F[free[i]];
    Solution S;
    CgResult cg;
    if (o.amg) {
        std::vector<Vec3> Xf;
        for (auto i : freenodes) Xf.push_back(M.X[i]);
        double h = o.aggregate_size;
        if (h <= 0) {   // three times the median edge of the first group
            std::vector<double> el;
            auto& G = M.groups.at(0);
            int k = nodes_per_element(G.kind);
            for (size_t e = 0; e < G.size(); e += std::max<size_t>(1, G.size() / 2000)) el.push_back(norm(M.X[G.conn[e * k + 1]] - M.X[G.conn[e * k]]));
            std::nth_element(el.begin(), el.begin() + el.size() / 2, el.end());
            h = 3.0 * el[el.size() / 2];
        }
        SaAmg amg;
        amg.setup(Kf, Xf, h, o.verbose);
        cg = pcg(Kf, Ff, uf, [&](const std::vector<double>& r, std::vector<double>& z) { amg.apply(r, z); }, o.rtol, o.maxiter);
    } else {
        auto d = Kf.diagonal();
        cg = pcg(Kf, Ff, uf, [&](const std::vector<double>& r, std::vector<double>& z) {
            z.resize(r.size());
            for (size_t i = 0; i < r.size(); i++) z[i] = r[i] / d[i];
        }, o.rtol, o.maxiter);
    }
    if (!cg.converged) throw std::runtime_error("PCG did not converge (residual " + std::to_string(cg.residual) + ")");
    S.iterations = cg.iterations; S.residual = cg.residual;
    std::vector<double> u(3 * nn, 0.0);
    for (size_t i = 0; i < free.size(); i++) u[free[i]] = uf[i];
    S.U.resize(nn);
    for (int64_t i = 0; i < nn; i++) S.U[i] = {u[3 * i], u[3 * i + 1], u[3 * i + 2]};
    std::vector<double> Ku;
    K.multiply(u, Ku);
    for (int64_t i = 0; i < nn; i++)
        if (isfix[i]) {
            Vec3 r{Ku[3 * i] - F[3 * i], Ku[3 * i + 1] - F[3 * i + 1], Ku[3 * i + 2] - F[3 * i + 2]};
            S.reaction += r;
            S.reaction_moment += cross(M.X[i], r);
        }
    // von Mises stress per element (largest at the Gauss points) and nodal stresses: the stress components
    // at the Gauss points are extrapolated to the corner nodes (trilinear for hexahedra, linear for
    // tetrahedra), mid-side nodes take the mean of their edge, and every node averages the elements
    // sharing it (the usual recovery, as in CalculiX); the nodal von Mises stress is that of the mean.
    double D[6][6];
    elasticity(E, nu, D);
    std::vector<std::array<double, 6>> acc(nn, std::array<double, 6>{});
    std::vector<double> cnt(nn, 0.0);
    double N[20], dn[20][3], dN[20][3];
    auto mises = [](const double* s) {
        return std::sqrt(0.5 * ((s[0] - s[1]) * (s[0] - s[1]) + (s[1] - s[2]) * (s[1] - s[2]) + (s[2] - s[0]) * (s[2] - s[0])) +
                         3 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]));
    };
    const int hex_edges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    const int tet_edges[6][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};
    for (auto& G : M.groups) {
        GaussRule g = rule(G.kind);
        int k = nodes_per_element(G.kind);
        bool tet = G.kind == ElementKind::Tet10;
        int nc = tet ? 4 : 8;
        // extrapolation weights W[c][q] from the Gauss points to the corner nodes
        std::vector<std::vector<double>> W(nc, std::vector<double>(g.p.size(), 0.0));
        if (tet) {
            double a = g.p[1][0], b = g.p[0][0];   // points near node q: volume coordinate a there, b elsewhere
            for (int c = 0; c < 4; c++)
                for (int q = 0; q < 4; q++) W[c][q] = ((c == q ? 1.0 : 0.0) - b) / (a - b);
        } else {
            for (int c = 0; c < 8; c++)
                for (size_t q = 0; q < g.p.size(); q++) {
                    if (g.p[q][0] == 0 || g.p[q][1] == 0 || g.p[q][2] == 0) continue;   // 27-point rule: corner points only
                    double w = 1;
                    for (int d = 0; d < 3; d++) w *= (1 + kHex20C[c][d] / std::fabs(g.p[q][d]) * (g.p[q][d] > 0 ? 1 : -1)) / 2;
                    W[c][q] = w;
                }
        }
        std::vector<double> vm(G.size(), 0.0);
        std::vector<std::array<double, 6>> sq(g.p.size());
        for (size_t e = 0; e < G.size(); e++) {
            const int64_t* conn = &G.conn[e * k];
            for (size_t q = 0; q < g.p.size(); q++) {
                shape_functions(G.kind, g.p[q][0], g.p[q][1], g.p[q][2], N, dn);
                gradients(M, conn, k, dn, dN);
                double eps[6] = {};
                for (int a = 0; a < k; a++) {
                    const Vec3& U = S.U[conn[a]];
                    eps[0] += dN[a][0] * U.x; eps[1] += dN[a][1] * U.y; eps[2] += dN[a][2] * U.z;
                    eps[3] += dN[a][1] * U.x + dN[a][0] * U.y; eps[4] += dN[a][2] * U.y + dN[a][1] * U.z;
                    eps[5] += dN[a][2] * U.x + dN[a][0] * U.z;
                }
                for (int r = 0; r < 6; r++) { sq[q][r] = 0; for (int l = 0; l < 6; l++) sq[q][r] += D[r][l] * eps[l]; }
                vm[e] = std::max(vm[e], mises(sq[q].data()));
            }
            double sn[20][6] = {};
            for (int c = 0; c < nc; c++)
                for (size_t q = 0; q < g.p.size(); q++)
                    if (W[c][q] != 0) for (int r = 0; r < 6; r++) sn[c][r] += W[c][q] * sq[q][r];
            if (k > nc) {
                const int (*ed)[2] = tet ? tet_edges : hex_edges;
                for (int m = 0; m < k - nc; m++)
                    for (int r = 0; r < 6; r++) sn[nc + m][r] = 0.5 * (sn[ed[m][0]][r] + sn[ed[m][1]][r]);
            }
            for (int a = 0; a < k; a++) { for (int r = 0; r < 6; r++) acc[conn[a]][r] += sn[a][r]; cnt[conn[a]] += 1; }
        }
        S.vm_elem.push_back(std::move(vm));
    }
    S.vm_node.resize(nn);
    for (int64_t i = 0; i < nn; i++) {
        double s[6];
        for (int r = 0; r < 6; r++) s[r] = acc[i][r] / std::max(cnt[i], 1.0);
        S.vm_node[i] = mises(s);
    }
    if (o.verbose)
        std::fprintf(stderr, "FE: %lld nodes, %zu elements, %d PCG iterations, residual %.2e, %.1f s\n", (long long)nn, M.n_elements(),
                     S.iterations, S.residual, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    return S;
}

// ------------------------------------------------------------------ meshing helpers
Mesh box_mesh(const Vec3& L, int nx, int ny, int nz, const Vec3& o, int order) {
    if (nx < 1 || ny < 1 || nz < 1) throw std::invalid_argument("box_mesh needs at least one element per direction");
    Mesh m;
    for (int k = 0; k <= nz; k++)
        for (int j = 0; j <= ny; j++)
            for (int i = 0; i <= nx; i++) m.X.push_back({o.x + L.x * i / nx, o.y + L.y * j / ny, o.z + L.z * k / nz});
    auto id = [&](int i, int j, int k) { return int64_t((k * (ny + 1) + j) * (nx + 1) + i); };
    ElementGroup g{ElementKind::Hex8, {}};
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                for (int64_t v : {id(i, j, k), id(i + 1, j, k), id(i + 1, j + 1, k), id(i, j + 1, k), id(i, j, k + 1), id(i + 1, j, k + 1),
                                  id(i + 1, j + 1, k + 1), id(i, j + 1, k + 1)})
                    g.conn.push_back(v);
    m.groups.push_back(g);
    return order == 2 ? to_hex20(m) : m;
}

Mesh to_hex20(const Mesh& m) {
    static const int edges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    Mesh o;
    o.X = m.X;
    std::map<std::pair<int64_t, int64_t>, int64_t> mid;
    for (auto& G : m.groups) {
        if (G.kind != ElementKind::Hex8) { o.groups.push_back(G); continue; }
        ElementGroup h{ElementKind::Hex20, {}};
        for (size_t e = 0; e < G.size(); e++) {
            const int64_t* c = &G.conn[e * 8];
            for (int a = 0; a < 8; a++) h.conn.push_back(c[a]);
            for (auto& ed : edges) {
                auto key = std::minmax(c[ed[0]], c[ed[1]]);
                auto it = mid.find(key);
                if (it == mid.end()) {
                    it = mid.emplace(key, int64_t(o.X.size())).first;
                    o.X.push_back((m.X[key.first] + m.X[key.second]) * 0.5);
                }
                h.conn.push_back(it->second);
            }
        }
        o.groups.push_back(h);
    }
    return o;
}

void write_vtk(const std::string& path, const Mesh& m, const Solution& s) {
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) throw std::runtime_error("cannot write " + path);
    std::fprintf(f, "# vtk DataFile Version 3.0\nTillBridge\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS %zu double\n", m.n_nodes());
    for (auto& x : m.X) std::fprintf(f, "%.9g %.9g %.9g\n", x.x, x.y, x.z);
    size_t ne = m.n_elements(), sz = 0;
    for (auto& G : m.groups) sz += G.size() * (nodes_per_element(G.kind) + 1);
    std::fprintf(f, "CELLS %zu %zu\n", ne, sz);
    // VTK node order equals the ANSYS/TillBridge order for hex8, hex20 and tet10
    for (auto& G : m.groups) {
        int k = nodes_per_element(G.kind);
        for (size_t e = 0; e < G.size(); e++) {
            std::fprintf(f, "%d", k);
            for (int a = 0; a < k; a++) std::fprintf(f, " %lld", (long long)G.conn[e * k + a]);
            std::fprintf(f, "\n");
        }
    }
    std::fprintf(f, "CELL_TYPES %zu\n", ne);
    for (auto& G : m.groups) {
        int t = G.kind == ElementKind::Hex8 ? 12 : (G.kind == ElementKind::Tet10 ? 24 : 25);
        for (size_t e = 0; e < G.size(); e++) std::fprintf(f, "%d\n", t);
    }
    std::fprintf(f, "POINT_DATA %zu\nVECTORS displacement double\n", m.n_nodes());
    for (auto& u : s.U) std::fprintf(f, "%.9g %.9g %.9g\n", u.x, u.y, u.z);
    std::fprintf(f, "SCALARS von_mises double 1\nLOOKUP_TABLE default\n");
    for (double v : s.vm_node) std::fprintf(f, "%.9g\n", v);
    if (!s.vm_elem.empty()) {   // element values: largest von Mises stress at the Gauss points
        std::fprintf(f, "CELL_DATA %zu\nSCALARS von_mises_element double 1\nLOOKUP_TABLE default\n", ne);
        for (const auto& g : s.vm_elem) for (double v : g) std::fprintf(f, "%.9g\n", v);
    }
    std::fclose(f);
}

void write_calculix(const std::string& path, const Mesh& m, double E, double nu, const std::vector<int64_t>& fixed,
                    const std::map<int64_t, Vec3>& loads) {
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) throw std::runtime_error("cannot write " + path);
    std::fprintf(f, "** Written by TillBridge: linear static model, units of the mesh (e.g. mm, N, MPa)\n*NODE, NSET=NALL\n");
    for (size_t i = 0; i < m.n_nodes(); i++) std::fprintf(f, "%zu, %.10g, %.10g, %.10g\n", i + 1, m.X[i].x, m.X[i].y, m.X[i].z);
    long id = 1;
    for (const auto& g : m.groups) {
        const char* type = g.kind == ElementKind::Tet10 ? "C3D10" : (g.kind == ElementKind::Hex8 ? "C3D8" : (g.kind == ElementKind::Hex20 ? "C3D20" : "C3D20R"));
        int nn = nodes_per_element(g.kind);
        std::fprintf(f, "*ELEMENT, TYPE=%s, ELSET=EALL\n", type);
        for (size_t e = 0; e < g.size(); e++) {
            std::fprintf(f, "%ld", id++);
            for (int k = 0; k < nn; k++)   // CalculiX lines hold at most 16 entries: continue after the 15th node
                std::fprintf(f, "%s%lld", k == 15 ? ",\n" : ", ", static_cast<long long>(g.conn[e * nn + k] + 1));
            std::fprintf(f, "\n");
        }
    }
    std::fprintf(f, "*NSET, NSET=FIXED\n");
    for (size_t k = 0; k < fixed.size(); k++) std::fprintf(f, "%lld%s", static_cast<long long>(fixed[k] + 1), (k % 10 == 9 || k + 1 == fixed.size()) ? "\n" : ", ");
    std::fprintf(f, "*NSET, NSET=LOADED\n");
    size_t k = 0;
    for (const auto& l : loads) { std::fprintf(f, "%lld%s", static_cast<long long>(l.first + 1), (k % 10 == 9 || k + 1 == loads.size()) ? "\n" : ", "); k++; }
    std::fprintf(f, "*MATERIAL, NAME=STEEL\n*ELASTIC\n%.10g, %.10g\n*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL\n", E, nu);
    std::fprintf(f, "*STEP\n*STATIC\n*BOUNDARY\nFIXED, 1, 3\n*CLOAD\n");
    for (const auto& l : loads)
        for (int d = 0; d < 3; d++)
            if (l.second[d] != 0) std::fprintf(f, "%lld, %d, %.12g\n", static_cast<long long>(l.first + 1), d + 1, l.second[d]);
    std::fprintf(f, "*NODE FILE\nU\n*EL FILE\nS\n*NODE PRINT, NSET=LOADED\nU\n*NODE PRINT, NSET=FIXED, TOTALS=ONLY\nRF\n*EL PRINT, ELSET=EALL\nS\n*END STEP\n");
    std::fclose(f);
}

}  // namespace tb
