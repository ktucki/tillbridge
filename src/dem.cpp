// Copyright (c) 2026 Karol Tucki. MIT License.
#include "tillbridge/dem.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <new>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace tb {

constexpr int kMaxContacts = 24;   // particle-particle contacts kept per sphere
constexpr Vec3 wn[5] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}};   // wall normals

// ------------------------------------------------------------------ contact law
double damping_beta(double e) {
    if (e >= 1.0) return 0.0;
    double l = std::log(e);
    return l / std::sqrt(l * l + kPi * kPi);
}
double effective_modulus(double E1, double nu1, double E2, double nu2) {
    return 1.0 / ((1 - nu1 * nu1) / E1 + (1 - nu2 * nu2) / E2);
}
double effective_shear_modulus(double E1, double nu1, double E2, double nu2) {
    double G1 = E1 / (2 * (1 + nu1)), G2 = E2 / (2 * (1 + nu2));
    return 1.0 / ((2 - nu1) / G1 + (2 - nu2) / G2);
}
double rayleigh_time(double R, double rho, double E, double nu) {
    double G = E / (2 * (1 + nu));
    return kPi * R * std::sqrt(rho / G) / (0.1631 * nu + 0.8766);
}

Vec3 hertz_mindlin(double del, const Vec3& n, const Vec3& vrel, double Reff, double meff, double Eeff, double Geff,
                   double mu, double beta, double dt, double k_coh, double area, Vec3& ut) {
    double sq = std::sqrt(Reff * del);
    double kn = 4.0 / 3.0 * Eeff * sq, Sn = 2.0 * Eeff * sq, St = 8.0 * Geff * sq;
    double vn = dot(vrel, n);
    Vec3 vt = vrel - n * vn;
    double fn = kn * del + 2.0 * std::sqrt(5.0 / 6.0) * beta * std::sqrt(Sn * meff) * vn;
    if (fn < 0) fn = 0;
    double fc = k_coh * area;
    ut = ut + vt * dt;
    ut = ut - n * dot(ut, n);
    Vec3 ft = ut * (-St) + vt * (2.0 * std::sqrt(5.0 / 6.0) * beta * std::sqrt(St * meff));
    double ftm = norm(ft), lim = mu * fn;
    if (ftm > lim && ftm > 0) { ft = ft * (lim / ftm); ut = ft * (-1.0 / St); }
    return n * (fn - fc) + ft;
}

Vec3 closest_on_triangle(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c) {
    Vec3 ab = b - a, ac = c - a, ap = p - a;
    double d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0 && d2 <= 0) return a;
    Vec3 bp = p - b;
    double d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) return b;
    double vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) return a + ab * (d1 / (d1 - d3));
    Vec3 cp = p - c;
    double d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) return c;
    double vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) return a + ac * (d2 / (d2 - d6));
    double va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
    double den = 1.0 / (va + vb + vc);
    return a + ab * (vb * den) + ac * (vc * den);
}

// ------------------------------------------------------------------ solver state
struct Contact {
    int j;
    Vec3 ut, mr;     // tangential spring, rolling moment
};
constexpr Contact kNoContact{-1, Vec3{}, Vec3{}};

struct Simulation::Impl {
    Params P;
    int NB = 0, NS = 0, NT = 0;
    // bodies
    std::vector<Vec3> bx, bv, bw, bf, bt;
    std::vector<Quat> bq;
    std::vector<double> bm, bI;
    std::vector<int> bfirst, bns, bmat;
    // spheres
    std::vector<Vec3> srel, sx, sf, stq;
    std::vector<double> sr, sE;
    std::vector<int> sb, snc, sntc;
    std::vector<Contact> sc;            // NS * kMaxContacts
    std::vector<Contact> stc;           // tool contact history
    std::vector<Vec3> swut;             // NS * 5 wall histories
    // tool
    std::vector<Vec3> ta, tbv, tc;
    // grids
    double cs = 0, ox = -0.05, oy = -0.05, oz = -0.05;
    int gx = 1, gy = 1, gz = 1;
    std::vector<int> head, nxt;
    double tox = 0, toy = 0, toz = 0, tcs = 0;
    int tgx = 1, tgy = 1, tgz = 1;
    std::vector<int> tstart, tlist;
    double rmax = 0;
    // bonds
    std::vector<int> bs_, bj_;
    std::vector<char> balive;
    std::vector<double> bfn, bmt;
    std::vector<Vec3> bfs, bmb;
    long nbroken_t = 0, nbroken_s = 0;
    // departures
    std::vector<Vec3> dv, dxp;
    std::vector<double> dtt;

    void setup(const Bed& bed, const Tool& tool);
    void build_tool_grid();
    void create_bonds();
    void bin_particles();
};

void Simulation::Impl::setup(const Bed& bed, const Tool& tool) {
    NB = int(bed.size()); NS = int(bed.n_spheres()); NT = int(tool.size());
    if (NB <= 0) throw std::invalid_argument("empty bed");
    bx.resize(NB); bv.resize(NB); bw.resize(NB); bf.resize(NB); bt.resize(NB); bq.assign(NB, Quat{});
    bm.resize(NB); bI.resize(NB); bfirst.resize(NB); bns.resize(NB); bmat.resize(NB);
    srel.resize(NS); sx.resize(NS); sf.resize(NS); stq.resize(NS); sr.resize(NS); sE.resize(NS); sb.resize(NS);
    snc.assign(NS, 0); sntc.assign(NS, 0); sc.assign(size_t(NS) * kMaxContacts, kNoContact); stc.assign(NS, kNoContact);
    swut.assign(size_t(NS) * 5, Vec3{});
    int k = 0;
    for (int b = 0; b < NB; b++) {
        const Body& B = bed.bodies[b];
        if (B.m <= 0 || B.I <= 0 || B.spheres.empty()) throw std::invalid_argument("bad body " + std::to_string(b));
        bx[b] = B.x; bv[b] = B.v; bw[b] = B.w; bm[b] = B.m; bI[b] = B.I; bmat[b] = B.mat;
        bfirst[b] = k; bns[b] = int(B.spheres.size());
        for (auto& s : B.spheres) {
            if (s.r <= 0) throw std::invalid_argument("bad sphere radius in body " + std::to_string(b));
            srel[k] = s.offset; sr[k] = s.r; sb[k] = b;
            sE[k] = (B.mat == kStone && P.E_stone > 0) ? P.E_stone : P.E;
            k++;
        }
    }
    ta.resize(NT); tbv.resize(NT); tc.resize(NT);
    for (int t = 0; t < NT; t++) { ta[t] = tool.tri[t][0]; tbv[t] = tool.tri[t][1]; tc[t] = tool.tri[t][2]; }
    for (int s = 0; s < NS; s++) rmax = std::max(rmax, sr[s]);
    cs = 2.05 * rmax;
    gx = int((P.Lx + 0.1) / cs) + 1; gy = int((P.Ly + 0.1) / cs) + 1;
    double zmax = P.Lz;
    for (int b = 0; b < NB; b++) if (bx[b].z + rmax * 2 > zmax) zmax = bx[b].z + rmax * 2;
    gz = int((zmax + 0.1) / cs) + 1;
    head.assign(size_t(gx) * gy * gz, -1); nxt.assign(NS, -1);
    dv.assign(NS, Vec3{}); dxp.assign(NS, Vec3{}); dtt.assign(NS, -1.0);
    build_tool_grid();
}

// Buckets of the tool frame (CSR): every bucket lists the triangles within rmax of it.
void Simulation::Impl::build_tool_grid() {
    double tmx = -1e30, tmy = -1e30, tmz = -1e30;
    tox = toy = toz = 1e30;
    for (int t = 0; t < NT; t++)
        for (const Vec3* p : {&ta[t], &tbv[t], &tc[t]}) {
            tox = std::fmin(tox, p->x); tmx = std::fmax(tmx, p->x); toy = std::fmin(toy, p->y);
            tmy = std::fmax(tmy, p->y); toz = std::fmin(toz, p->z); tmz = std::fmax(tmz, p->z);
        }
    tcs = std::fmax(2.0 * rmax, 0.02);
    tgx = tgy = tgz = 1;
    if (NT > 0) {
        tox -= tcs; toy -= tcs; toz -= tcs;
        tgx = int((tmx - tox) / tcs) + 3; tgy = int((tmy - toy) / tcs) + 3; tgz = int((tmz - toz) / tcs) + 3;
    }
    int ncell = tgx * tgy * tgz;
    tstart.assign(ncell + 1, 0);
    std::vector<int> fill(ncell, 0);
    for (int pass = 0; pass < 2; pass++) {
        for (int t = 0; t < NT; t++) {
            double x0 = std::fmin(ta[t].x, std::fmin(tbv[t].x, tc[t].x)) - rmax, x1 = std::fmax(ta[t].x, std::fmax(tbv[t].x, tc[t].x)) + rmax;
            double y0 = std::fmin(ta[t].y, std::fmin(tbv[t].y, tc[t].y)) - rmax, y1 = std::fmax(ta[t].y, std::fmax(tbv[t].y, tc[t].y)) + rmax;
            double z0 = std::fmin(ta[t].z, std::fmin(tbv[t].z, tc[t].z)) - rmax, z1 = std::fmax(ta[t].z, std::fmax(tbv[t].z, tc[t].z)) + rmax;
            int i0 = int((x0 - tox) / tcs), i1 = int((x1 - tox) / tcs), j0 = int((y0 - toy) / tcs), j1 = int((y1 - toy) / tcs);
            int k0 = int((z0 - toz) / tcs), k1 = int((z1 - toz) / tcs);
            for (int a = i0; a <= i1; a++)
                for (int bb = j0; bb <= j1; bb++)
                    for (int c = k0; c <= k1; c++) {
                        if (a < 0 || bb < 0 || c < 0 || a >= tgx || bb >= tgy || c >= tgz) continue;
                        int id = (c * tgy + bb) * tgx + a;
                        if (pass == 0) tstart[id + 1]++;
                        else tlist[tstart[id] + fill[id]++] = t;
                    }
        }
        if (pass == 0) {
            for (int i = 0; i < ncell; i++) tstart[i + 1] += tstart[i];
            tlist.assign(tstart[ncell] + 1, 0);
        }
    }
}

void Simulation::Impl::bin_particles() {
    std::fill(head.begin(), head.end(), -1);
    for (int s = 0; s < NS; s++) {
        int ix = int((sx[s].x - ox) / cs), iy = int((sx[s].y - oy) / cs), iz = int((sx[s].z - oz) / cs);
        ix = ix < 0 ? 0 : (ix >= gx ? gx - 1 : ix);
        iy = iy < 0 ? 0 : (iy >= gy ? gy - 1 : iy);
        iz = iz < 0 ? 0 : (iz >= gz ? gz - 1 : iz);
        int id = (iz * gy + iy) * gx + ix;
        nxt[s] = head[id]; head[id] = s;
    }
}

// Bonds between spheres of different bodies closer than (r_i + r_j)(1 + bond_gap), at t = 0.
void Simulation::Impl::create_bonds() {
    for (int s = 0; s < NS; s++) {
        int a0 = int((sx[s].x - ox) / cs), b0 = int((sx[s].y - oy) / cs), c0 = int((sx[s].z - oz) / cs);
        for (int dz = -1; dz <= 1; dz++)
            for (int dy = -1; dy <= 1; dy++)
                for (int dx = -1; dx <= 1; dx++) {
                    int a = a0 + dx, bb = b0 + dy, c = c0 + dz;
                    if (a < 0 || bb < 0 || c < 0 || a >= gx || bb >= gy || c >= gz) continue;
                    for (int j = head[(c * gy + bb) * gx + a]; j >= 0; j = nxt[j]) {
                        if (j <= s || sb[j] == sb[s]) continue;
                        if (bmat[sb[s]] == kFixed && bmat[sb[j]] == kFixed) continue;
                        if (bx[sb[s]].z < P.bond_zmin || bx[sb[j]].z < P.bond_zmin) continue;
                        double dist = norm(sx[s] - sx[j]);
                        if (dist > (sr[s] + sr[j]) * (1.0 + P.bond_gap)) continue;
                        bs_.push_back(s); bj_.push_back(j);
                    }
                }
    }
    size_t n = bs_.size();
    balive.assign(n, 1); bfn.assign(n, 0); bmt.assign(n, 0); bfs.assign(n, Vec3{}); bmb.assign(n, Vec3{});
    std::fprintf(stderr, "tillbridge: %zu parallel bonds created\n", n);
}

Simulation::Simulation(const Params& p, const Bed& bed, const Tool& tool) : d_(new Impl) {
    p.check();
    d_->P = p;
    d_->setup(bed, tool);
}
Simulation::~Simulation() { delete d_; }

RunResult Simulation::run(bool verbose) {
    Impl& S = *d_;
    const Params& P = S.P;
    const int NB = S.NB, NS = S.NS, NT = S.NT;
    RunResult R;
    const double beta = damping_beta(P.e);
    const Vec3 toolv{P.v_tool, 0, 0};
    std::vector<double> tF(5 * size_t(NT + 1), 0.0), tFacc(5 * size_t(NT + 1), 0.0);
    std::vector<double> tW(2 * size_t(NT + 1), 0.0), tWacc(2 * size_t(NT + 1), 0.0);
    long nacc = 0;
    double Facc[3] = {0, 0, 0};
    long long ncontact_pp = 0;
    auto c0 = std::chrono::steady_clock::now();
    double t = 0, nextout = 0;
    long step = 0;
    double kemin = 1e300, mtot = 0;
    // energy accounting over the averaging window (Section "energy closure" of docs/model.md):
    // work of the tool, energy taken up by the tool contacts (relative motion), by all other contacts
    // (particle-particle, walls, rolling resistance, bonds), and the change of kinetic and potential energy
    double eTool = 0, eToolRel = 0, eBulk = 0, eMech0 = 0, eMech1 = 0, ePot0 = 0, ePot1 = 0;
    bool winStarted = false;
    auto mech = [&](double& kin, double& pot) {
        kin = 0; pot = 0;
        for (int b = 0; b < NB; b++) {
            kin += 0.5 * S.bm[b] * dot(S.bv[b], S.bv[b]) + 0.5 * S.bI[b] * dot(S.bw[b], S.bw[b]);
            pot += S.bm[b] * P.g * S.bx[b].z;
        }
    };
    for (int b = 0; b < NB; b++) mtot += S.bm[b];
    // effective moduli are constant per pair of materials
    const double Eww = effective_modulus(P.E, P.nu, P.E_wall, P.nu_wall), Gww = effective_shear_modulus(P.E, P.nu, P.E_wall, P.nu_wall);

    while (t < P.t_end - 0.5 * P.dt) {
        const bool moving = t >= P.t_tool;
        const Vec3 tp = moving ? toolv * (t - P.t_tool) : Vec3{};
        const bool inWin = moving && tp.x >= P.avg_x0 && tp.x <= P.avg_x1;
        if (inWin && !winStarted) { mech(eMech0, ePot0); winStarted = true; }
#pragma omp parallel for
        for (int s = 0; s < NS; s++) { int b = S.sb[s]; S.sx[s] = S.bx[b] + rotate(S.bq[b], S.srel[s]); }
        S.bin_particles();
        if (step == 0 && P.bond_E > 0) S.create_bonds();
        std::fill(tF.begin(), tF.end(), 0.0);
        std::fill(tW.begin(), tW.end(), 0.0);
        double ftx = 0, fty = 0, ftz = 0;
        long long ncp = 0;
        double pBulk = 0, pTool = 0, pToolRel = 0;
#pragma omp parallel for schedule(dynamic, 256) reduction(+ : ftx, fty, ftz, ncp, pBulk, pTool, pToolRel)
        for (int s = 0; s < NS; s++) {
            Vec3 F, T;
            const int b = S.sb[s];
            alignas(Contact) unsigned char newc_buf[sizeof(Contact) * kMaxContacts];
            Contact* newc = reinterpret_cast<Contact*>(newc_buf);
            int nn = 0;
            Contact* old = &S.sc[size_t(s) * kMaxContacts];
            const int a0 = int((S.sx[s].x - S.ox) / S.cs), b0 = int((S.sx[s].y - S.oy) / S.cs), cz0 = int((S.sx[s].z - S.oz) / S.cs);
            for (int dz = -1; dz <= 1; dz++)
                for (int dy = -1; dy <= 1; dy++)
                    for (int dx = -1; dx <= 1; dx++) {
                        int a = a0 + dx, bb = b0 + dy, c = cz0 + dz;
                        if (a < 0 || bb < 0 || c < 0 || a >= S.gx || bb >= S.gy || c >= S.gz) continue;
                        for (int j = S.head[(c * S.gy + bb) * S.gx + a]; j >= 0; j = S.nxt[j]) {
                            if (S.sb[j] == b) continue;
                            Vec3 d = S.sx[s] - S.sx[j];
                            double dist = norm(d), del = S.sr[s] + S.sr[j] - dist;
                            if (del <= 0 || dist <= 0) continue;
                            Vec3 n = d * (1.0 / dist);
                            int bj = S.sb[j];
                            // velocities of both bodies at the contact point (middle of the overlap)
                            Vec3 xc = S.sx[s] - n * (S.sr[s] - 0.5 * del);
                            Vec3 vi = S.bv[b] + cross(S.bw[b], xc - S.bx[b]);
                            Vec3 vj = S.bv[bj] + cross(S.bw[bj], xc - S.bx[bj]);
                            Vec3 ut, mr;
                            for (int q = 0; q < S.snc[s]; q++) if (old[q].j == j) { ut = old[q].ut; mr = old[q].mr; break; }
                            double Reff = S.sr[s] * S.sr[j] / (S.sr[s] + S.sr[j]), meff = S.bm[b] * S.bm[bj] / (S.bm[b] + S.bm[bj]);
                            double Ee = effective_modulus(S.sE[s], P.nu, S.sE[j], P.nu);
                            Vec3 Fc = hertz_mindlin(del, n, vi - vj, Reff, meff, Ee, effective_shear_modulus(S.sE[s], P.nu, S.sE[j], P.nu),
                                                    P.mu_pp, beta, P.dt, P.k_coh, kPi * Reff * del, ut);
                            ncp++;
                            F = F + Fc;
                            T = T + cross(xc - S.bx[b], Fc);
                            // rolling resistance: k_r = 2.25 k_n mu_r^2 R*^2 with the tangent normal stiffness
                            // k_n = 2 E* (R* del)^1/2, moment limited to mu_r R* F_n; twisting excluded
                            if (P.mu_r > 0) {
                                double knt = 2.0 * Ee * std::sqrt(Reff * del);
                                double fnn = dot(Fc, n) + P.k_coh * kPi * Reff * del;   // repulsive part
                                Vec3 wr = S.bw[b] - S.bw[bj];
                                wr = wr - n * dot(wr, n);
                                mr = mr - n * dot(mr, n);
                                double kr = 2.25 * knt * P.mu_r * P.mu_r * Reff * Reff;
                                mr = mr - wr * (kr * P.dt);
                                double mlim = P.mu_r * Reff * (fnn > 0 ? fnn : 0), mm = norm(mr);
                                if (mm > mlim) mr = mm > 0 ? mr * (mlim / mm) : mr;
                                T = T + mr;
                                pBulk -= dot(mr, S.bw[b]);
                            }
                            pBulk -= dot(Fc, vi);
                            if (nn < kMaxContacts) { new (&newc[nn]) Contact{j, ut, mr}; nn++; }
                        }
                    }
            std::copy(newc, newc + nn, old);
            S.snc[s] = nn;
            // box walls x = 0, x = Lx, y = 0, y = Ly and floor z = 0, with tangential history
            const double wd[5] = {S.sx[s].x, P.Lx - S.sx[s].x, S.sx[s].y, P.Ly - S.sx[s].y, S.sx[s].z};
            for (int w = 0; w < 5; w++) {
                Vec3& hist = S.swut[size_t(s) * 5 + w];
                double del = S.sr[s] - wd[w];
                if (del <= 0) { hist = Vec3{}; continue; }
                Vec3 ut = hist;
                Vec3 xc = S.sx[s] - wn[w] * (S.sr[s] - del);
                Vec3 vi = S.bv[b] + cross(S.bw[b], xc - S.bx[b]);
                double Ew = S.sE[s] == P.E ? Eww : effective_modulus(S.sE[s], P.nu, P.E_wall, P.nu_wall);
                double Gw = S.sE[s] == P.E ? Gww : effective_shear_modulus(S.sE[s], P.nu, P.E_wall, P.nu_wall);
                Vec3 Fc = hertz_mindlin(del, wn[w], vi, S.sr[s], S.bm[b], Ew, Gw, P.mu_pw, beta, P.dt, 0.0, 0.0, ut);
                hist = ut;
                pBulk -= dot(Fc, vi);
                F = F + Fc;
                T = T + cross(xc - S.bx[b], Fc);
            }
            // tool: the nearest triangle of the bucket
            S.sntc[s] = 0;
            if (moving && NT > 0) {
                Vec3 pl = S.sx[s] - tp;
                int ia = int(std::floor((pl.x - S.tox) / S.tcs)), ib = int(std::floor((pl.y - S.toy) / S.tcs)),
                    ic = int(std::floor((pl.z - S.toz) / S.tcs));
                if (ia >= 0 && ib >= 0 && ic >= 0 && ia < S.tgx && ib < S.tgy && ic < S.tgz) {
                    int id = (ic * S.tgy + ib) * S.tgx + ia;
                    double best = S.sr[s];
                    int bt = -1;
                    Vec3 bcp;
                    for (int q = S.tstart[id]; q < S.tstart[id + 1]; q++) {
                        int tt = S.tlist[q];
                        Vec3 cp = closest_on_triangle(pl, S.ta[tt], S.tbv[tt], S.tc[tt]);
                        double dd = norm(pl - cp);
                        if (dd < best) { best = dd; bt = tt; bcp = cp; }
                    }
                    if (bt >= 0 && best > 0) {
                        double del = S.sr[s] - best;
                        Vec3 n = (pl - bcp) * (1.0 / best);
                        Vec3 ut;
                        if (S.stc[s].j >= 0) ut = S.stc[s].ut;   // history kept when the contact moves to a neighbouring facet
                        Vec3 xc = bcp + tp;
                        Vec3 vi = S.bv[b] + cross(S.bw[b], xc - S.bx[b]);
                        double Ew = S.sE[s] == P.E ? Eww : effective_modulus(S.sE[s], P.nu, P.E_wall, P.nu_wall);
                        double Gw = S.sE[s] == P.E ? Gww : effective_shear_modulus(S.sE[s], P.nu, P.E_wall, P.nu_wall);
                        Vec3 Fc = hertz_mindlin(del, n, vi - toolv, S.sr[s], S.bm[b], Ew, Gw, P.mu_pw, beta, P.dt, 0.0, 0.0, ut);
                        S.stc[s].j = bt; S.stc[s].ut = ut; S.sntc[s] = 1;
                        F = F + Fc;
                        T = T + cross(xc - S.bx[b], Fc);
                        ftx -= Fc.x; fty -= Fc.y; ftz -= Fc.z;
                        pTool += dot(Fc, toolv);
                        pToolRel -= dot(Fc, vi - toolv);
                        double fnrm = dot(Fc, n);
                        // dissipated power at the contact: tangential (friction) and normal (damping) part
                        Vec3 vr = vi - toolv;
                        double vn = dot(vr, n);
                        Vec3 vt = vr - n * vn, Ft = Fc - n * fnrm;
                        double pt = -dot(Ft, vt), pn = -fnrm * vn;
                        double* f5 = &tF[5 * size_t(bt)];
                        double* w2 = &tW[2 * size_t(bt)];
#pragma omp atomic
                        f5[0] += fnrm;
#pragma omp atomic
                        f5[1] += 1.0;
#pragma omp atomic
                        f5[2] -= Fc.x;
#pragma omp atomic
                        f5[3] -= Fc.y;
#pragma omp atomic
                        f5[4] -= Fc.z;
#pragma omp atomic
                        w2[0] += pt;
#pragma omp atomic
                        w2[1] += pn;
                    }
                }
            }
            if (!S.sntc[s] && S.stc[s].j >= 0) { S.dv[s] = S.bv[b]; S.dxp[s] = S.sx[s]; S.dtt[s] = tp.x; }
            if (!S.sntc[s]) S.stc[s].j = -1;
            S.sf[s] = F; S.stq[s] = T;
        }
        ncontact_pp += ncp;
        std::fill(S.bf.begin(), S.bf.end(), Vec3{});
        std::fill(S.bt.begin(), S.bt.end(), Vec3{});
        for (int s = 0; s < NS; s++) { int b = S.sb[s]; S.bf[b] = S.bf[b] + S.sf[s]; S.bt[b] = S.bt[b] + S.stq[s]; }
        // parallel bonds
        const int nbond = int(S.bs_.size());
        if (nbond) {
            long nbt = 0, nbs = 0;
            double pBond = 0;
#pragma omp parallel for schedule(static) reduction(+ : nbt, nbs, pBond)
            for (int q = 0; q < nbond; q++) {
                if (!S.balive[q]) continue;
                int s = S.bs_[q], j = S.bj_[q], bi = S.sb[s], bjj = S.sb[j];
                Vec3 d = S.sx[s] - S.sx[j];
                double dist = norm(d);
                if (dist <= 0) continue;
                Vec3 n = d * (1.0 / dist);
                double Rb = P.bond_lam * std::fmin(S.sr[s], S.sr[j]), A = kPi * Rb * Rb, I = 0.25 * kPi * Rb * Rb * Rb * Rb, J = 2 * I;
                double L = S.sr[s] + S.sr[j], kn = P.bond_E / L, ks = kn / P.bond_kr;
                Vec3 xc = S.sx[j] + n * (dist * S.sr[j] / L);
                Vec3 vi = S.bv[bi] + cross(S.bw[bi], xc - S.bx[bi]), vj = S.bv[bjj] + cross(S.bw[bjj], xc - S.bx[bjj]);
                Vec3 vr = vi - vj;
                double vn = dot(vr, n);
                Vec3 vt = vr - n * vn;
                Vec3 Fs = S.bfs[q]; Fs = Fs - n * dot(Fs, n);
                Vec3 Mb = S.bmb[q]; Mb = Mb - n * dot(Mb, n);
                S.bfn[q] -= kn * A * vn * P.dt;
                Fs = Fs - vt * (ks * A * P.dt);
                Vec3 wr = S.bw[bi] - S.bw[bjj];
                double wn = dot(wr, n);
                Vec3 wb = wr - n * wn;
                S.bmt[q] -= ks * J * wn * P.dt;
                Mb = Mb - wb * (kn * I * P.dt);
                S.bfs[q] = Fs; S.bmb[q] = Mb;
                double sig = -S.bfn[q] / A + norm(Mb) * Rb / I, tau = norm(Fs) / A + std::fabs(S.bmt[q]) * Rb / J;
                if (sig > P.bond_sc) { S.balive[q] = 0; nbt++; continue; }
                if (tau > P.bond_tc) { S.balive[q] = 0; nbs++; continue; }
                double me = S.bm[bi] * S.bm[bjj] / (S.bm[bi] + S.bm[bjj]);
                if (S.bmat[bi] == kFixed) me = S.bm[bjj];
                if (S.bmat[bjj] == kFixed) me = S.bm[bi];
                double cn = 2 * P.bond_zeta * std::sqrt(kn * A * me), ct = 2 * P.bond_zeta * std::sqrt(ks * A * me);
                Vec3 Fb = n * (S.bfn[q] - cn * vn) + Fs + vt * (-ct);
                Vec3 Ti = cross(xc - S.bx[bi], Fb) + Mb + n * S.bmt[q];
                Vec3 Tj = cross(xc - S.bx[bjj], Fb * -1.0) - (Mb + n * S.bmt[q]);
                pBond -= dot(Fb, S.bv[bi]) + dot(Ti, S.bw[bi]) - dot(Fb, S.bv[bjj]) + dot(Tj, S.bw[bjj]);
#pragma omp atomic
                S.bf[bi].x += Fb.x;
#pragma omp atomic
                S.bf[bi].y += Fb.y;
#pragma omp atomic
                S.bf[bi].z += Fb.z;
#pragma omp atomic
                S.bt[bi].x += Ti.x;
#pragma omp atomic
                S.bt[bi].y += Ti.y;
#pragma omp atomic
                S.bt[bi].z += Ti.z;
#pragma omp atomic
                S.bf[bjj].x -= Fb.x;
#pragma omp atomic
                S.bf[bjj].y -= Fb.y;
#pragma omp atomic
                S.bf[bjj].z -= Fb.z;
#pragma omp atomic
                S.bt[bjj].x += Tj.x;
#pragma omp atomic
                S.bt[bjj].y += Tj.y;
#pragma omp atomic
                S.bt[bjj].z += Tj.z;
            }
            S.nbroken_t += nbt; S.nbroken_s += nbs;
            pBulk += pBond;
        }
        // semi-implicit Euler, spherical inertia, quaternion update
#pragma omp parallel for
        for (int b = 0; b < NB; b++) {
            if (S.bmat[b] == kFixed) { S.bv[b] = Vec3{}; S.bw[b] = Vec3{}; continue; }   // anchored subsoil
            Vec3 a = S.bf[b] * (1.0 / S.bm[b]);
            a.z -= P.g;
            S.bv[b] = S.bv[b] + a * P.dt;
            S.bx[b] = S.bx[b] + S.bv[b] * P.dt;
            S.bw[b] = S.bw[b] + S.bt[b] * (P.dt / S.bI[b]);
            double wm = norm(S.bw[b]);
            if (wm > 0) {
                double h = wm * P.dt / 2;
                Vec3 ax = S.bw[b] * (std::sin(h) / wm);
                Quat dq{std::cos(h), ax.x, ax.y, ax.z};
                S.bq[b] = dq * S.bq[b];
                normalize(S.bq[b]);
            }
        }
        if (inWin) {
            eTool += pTool * P.dt; eToolRel += pToolRel * P.dt; eBulk += pBulk * P.dt;
            mech(eMech1, ePot1);
            for (int q = 0; q < 5 * NT; q++) tFacc[q] += tF[q];
            for (int q = 0; q < 2 * NT; q++) tWacc[q] += tW[q];
            Facc[0] += ftx; Facc[1] += fty; Facc[2] += ftz;
            nacc++;
        }
        if (t >= nextout - 1e-12) {
            double ke = 0, vmx = 0;
            for (int b = 0; b < NB; b++) {
                ke += 0.5 * S.bm[b] * dot(S.bv[b], S.bv[b]) + 0.5 * S.bI[b] * dot(S.bw[b], S.bw[b]);
                vmx = std::max(vmx, norm(S.bv[b]));
            }
            if (P.snap_dt > 0 && std::fmod(t + 1e-9, P.snap_dt) < P.out_dt - 1e-9) {
                char sn[1200];
                std::snprintf(sn, sizeof sn, "%s_snap_%04d.txt", P.out.c_str(), int(std::floor(t / P.snap_dt + 0.5)));
                if (std::FILE* fs = std::fopen(sn, "w")) {
                    for (int b = 0; b < NB; b++)
                        std::fprintf(fs, "%.4f %.4f %.4f %.3f %.3f %.3f\n", S.bx[b].x, S.bx[b].y, S.bx[b].z, S.bv[b].x, S.bv[b].y, S.bv[b].z);
                    std::fclose(fs);
                }
            }
            // a resting bed must not gain energy before the tool moves (optional check)
            bool grow = false;
            if (P.ke_growth > 0 && !moving && t >= 0.02 - 1e-12) {
                kemin = std::min(kemin, ke);
                if (ke > P.ke_growth * kemin + 0.005 * mtot) grow = true;
            }
            if (!(vmx <= P.v_max) || grow) {   // also catches NaN
                R.unstable = true; R.t_unstable = t; R.v_max_found = vmx; R.ke_unstable = ke; R.steps = step;
                std::fprintf(stderr, "tillbridge: unstable run at t = %.4f s (max particle speed %g m/s, kinetic energy %g J); "
                                     "reduce dt or the cohesion\n", t, vmx, ke);
                return R;
            }
            R.history.push_back({t, ftx, fty, ftz, ke, tp.x});
            if (verbose && R.history.size() % 100 == 1)
                std::fprintf(stderr, "t = %.3f s  Fx = %.1f N  Ekin = %.3g J\n", t, ftx, ke);
            nextout += P.out_dt;
        }
        t += P.dt;
        step++;
    }
    R.cpu_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - c0).count();
    R.steps = step;
    R.n_avg = nacc;
    R.mean_contacts_per_step = step ? double(ncontact_pp) / step : 0.0;
    R.bonds_broken_tension = S.nbroken_t; R.bonds_broken_shear = S.nbroken_s;
    if (nacc) R.F_mean = {Facc[0] / nacc, Facc[1] / nacc, Facc[2] / nacc};
    if (nacc) {
        R.energy.tool = eTool; R.energy.tool_contacts = eToolRel; R.energy.bulk = eBulk;
        R.energy.kinetic = eMech1 - eMech0; R.energy.potential = ePot1 - ePot0;
        R.energy.distance = nacc * P.dt * P.v_tool;
    }
    R.triangles.resize(NT);
    for (int q = 0; q < NT; q++) {
        TriangleLoad& L = R.triangles[q];
        if (!nacc) continue;
        L.Fn = tFacc[5 * q] / nacc; L.contacts = tFacc[5 * q + 1] / nacc;
        L.F = {tFacc[5 * q + 2] / nacc, tFacc[5 * q + 3] / nacc, tFacc[5 * q + 4] / nacc};
        L.Pt = tWacc[2 * q] / nacc; L.Pn = tWacc[2 * q + 1] / nacc;
    }
    for (int s = 0; s < NS; s++)
        if (S.dtt[s] >= 0) R.departures.push_back({s, S.sb[s], S.dxp[s], S.dv[s], S.dtt[s]});
    for (size_t q = 0; q < S.bs_.size(); q++) R.bonds.push_back({S.sb[S.bs_[q]], S.sb[S.bj_[q]], S.balive[q] != 0});
    R.state.bodies.resize(NB);
    for (int b = 0; b < NB; b++) {
        Body& B = R.state.bodies[b];
        B.x = S.bx[b]; B.m = S.bm[b]; B.I = S.bI[b]; B.mat = S.bmat[b]; B.v = S.bv[b]; B.w = S.bw[b];
        for (int s = S.bfirst[b]; s < S.bfirst[b] + S.bns[b]; s++) B.spheres.push_back({rotate(S.bq[b], S.srel[s]), S.sr[s]});
    }
    return R;
}

// ------------------------------------------------------------------ output
static std::FILE* open_w(const std::string& fn) {
    std::FILE* f = std::fopen(fn.c_str(), "w");
    if (!f) throw std::runtime_error("cannot write " + fn);
    return f;
}

std::map<std::string, double> RunResult::summary() const {
    std::map<std::string, double> m;
    if (unstable) {
        m["unstable"] = 1; m["t_unstable"] = t_unstable; m["v_max_found"] = v_max_found; m["ke_unstable"] = ke_unstable;
        return m;
    }
    m["bodies"] = double(state.size()); m["spheres"] = double(state.n_spheres()); m["triangles"] = double(triangles.size());
    m["steps"] = double(steps); m["cpu_s"] = cpu_s; m["mean_contacts_per_step"] = mean_contacts_per_step;
    m["bonds"] = double(bonds.size()); m["bonds_broken_tension"] = double(bonds_broken_tension);
    m["bonds_broken_shear"] = double(bonds_broken_shear); m["n_avg"] = double(n_avg);
    m["Fx_mean"] = F_mean.x; m["Fy_mean"] = F_mean.y; m["Fz_mean"] = F_mean.z;
    m["E_tool"] = energy.tool; m["E_tool_contacts"] = energy.tool_contacts; m["E_bulk"] = energy.bulk;
    m["dE_kinetic"] = energy.kinetic; m["dE_potential"] = energy.potential; m["E_distance"] = energy.distance;
    m["E_closure"] = energy.closure();
    return m;
}

void RunResult::write(const std::string& prefix) const {
    std::FILE* f = open_w(prefix + "_force.csv");
    std::fprintf(f, "t,Fx,Fy,Fz,Ekin,x_tool\n");
    for (auto& h : history) std::fprintf(f, "%.6f,%.4f,%.4f,%.4f,%.6e,%.5f\n", h.t, h.Fx, h.Fy, h.Fz, h.Ekin, h.x_tool);
    std::fclose(f);
    if (unstable) {
        f = open_w(prefix + "_summary.txt");
        std::fprintf(f, "version = %s\nunstable = 1\nt_unstable = %.6f\nv_max_found = %g\nke_unstable = %g\n", kVersion, t_unstable,
                     v_max_found, ke_unstable);
        std::fclose(f);
        return;
    }
    f = open_w(prefix + "_state.txt");
    std::fprintf(f, "%zu %zu\n", state.size(), state.n_spheres());
    for (auto& B : state.bodies) {
        std::fprintf(f, "%.6f %.6f %.6f %.6e %.6e %zu %d\n", B.x.x, B.x.y, B.x.z, B.m, B.I, B.spheres.size(), B.mat);
        for (auto& s : B.spheres) std::fprintf(f, "%.6f %.6f %.6f %.6f\n", s.offset.x, s.offset.y, s.offset.z, s.r);
    }
    std::fclose(f);
    f = open_w(prefix + "_vel.txt");
    for (auto& B : state.bodies) std::fprintf(f, "%.9e %.9e %.9e %.9e %.9e %.9e\n", B.v.x, B.v.y, B.v.z, B.w.x, B.w.y, B.w.z);
    std::fclose(f);
    f = open_w(prefix + "_depart.txt");
    for (auto& d : departures)
        std::fprintf(f, "%ld %ld %.5f %.5f %.5f %.4f %.4f %.4f %.4f\n", d.sphere, d.body, d.x.x, d.x.y, d.x.z, d.v.x, d.v.y, d.v.z, d.x_tool);
    std::fclose(f);
    if (!bonds.empty()) {
        f = open_w(prefix + "_bonds.txt");
        for (auto& b : bonds) std::fprintf(f, "%ld %ld %d\n", b.body_i, b.body_j, int(b.alive));
        std::fclose(f);
    }
    f = open_w(prefix + "_tri.csv");
    std::fprintf(f, "tri,Fn_mean,contacts_mean,Fx_mean,Fy_mean,Fz_mean\n");
    for (size_t q = 0; q < triangles.size(); q++) {
        auto& L = triangles[q];
        std::fprintf(f, "%zu,%.6e,%.6e,%.6e,%.6e,%.6e\n", q, L.Fn, L.contacts, L.F.x, L.F.y, L.F.z);
    }
    std::fclose(f);
    f = open_w(prefix + "_triW.csv");
    std::fprintf(f, "tri,Pt_mean,Pn_mean\n");
    for (size_t q = 0; q < triangles.size(); q++) std::fprintf(f, "%zu,%.6e,%.6e\n", q, triangles[q].Pt, triangles[q].Pn);
    std::fclose(f);
    f = open_w(prefix + "_summary.txt");
    std::fprintf(f, "version = %s\nbodies = %zu\nspheres = %zu\ntriangles = %zu\nsteps = %ld\ncpu_s = %.2f\n", kVersion, state.size(),
                 state.n_spheres(), triangles.size(), steps, cpu_s);
    std::fprintf(f, "mean_contacts_per_step = %.1f\n", mean_contacts_per_step);
    std::fprintf(f, "bonds = %zu\nbonds_broken_tension = %ld\nbonds_broken_shear = %ld\n", bonds.size(), bonds_broken_tension,
                 bonds_broken_shear);
    std::fprintf(f, "n_avg = %ld\nFx_mean = %.6f\nFy_mean = %.6f\nFz_mean = %.6f\n", n_avg, F_mean.x, F_mean.y, F_mean.z);
    std::fprintf(f, "E_distance = %.6f\nE_tool = %.6e\nE_tool_contacts = %.6e\nE_bulk = %.6e\ndE_kinetic = %.6e\ndE_potential = %.6e\nE_closure = %.6e\n",
                 energy.distance, energy.tool, energy.tool_contacts, energy.bulk, energy.kinetic, energy.potential, energy.closure());
    std::fclose(f);
}

RunResult run_dem(const Params& p, const Bed& bed, const Tool& tool, bool write, bool verbose) {
    Simulation sim(p, bed, tool);
    RunResult r = sim.run(verbose);
    if (write) r.write(p.out);
    return r;
}

}  // namespace tb
