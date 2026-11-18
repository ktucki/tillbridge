// Copyright (c) 2026 Karol Tucki. MIT License.
#include "tillbridge/calibrate.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include "tillbridge/bed.hpp"

namespace tb {

// ------------------------------------------------------------------ statistics
double normal_pdf(double z) { return std::exp(-0.5 * z * z) / std::sqrt(2 * 3.14159265358979323846); }
double normal_cdf(double z) { return 0.5 * std::erfc(-z / std::sqrt(2.0)); }
double expected_improvement(double mu, double sd, double best) {
    sd = std::max(sd, 1e-12);
    double z = (best - mu) / sd;
    return (best - mu) * normal_cdf(z) + sd * normal_pdf(z);
}

std::vector<std::vector<double>> latin_hypercube(int n, int d, uint64_t seed) {
    Rng rng(seed + 7);
    std::vector<std::vector<double>> U(n, std::vector<double>(d));
    for (int k = 0; k < d; k++) {
        std::vector<int> perm(n);
        std::iota(perm.begin(), perm.end(), 0);
        for (int i = n - 1; i > 0; i--) std::swap(perm[i], perm[rng.integer(0, i + 1)]);
        for (int i = 0; i < n; i++) U[i][k] = (perm[i] + rng.uniform()) / n;
    }
    return U;
}

std::vector<std::vector<double>> halton(int n, int d, uint64_t seed) {
    static const int primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71};
    if (d > 20) throw std::invalid_argument("at most 20 parameters");
    Rng rng(seed + 11);
    std::vector<double> shift(d);
    for (auto& s : shift) s = rng.uniform();
    std::vector<std::vector<double>> P(n, std::vector<double>(d));
    for (int i = 0; i < n; i++)
        for (int k = 0; k < d; k++) {
            double f = 1, r = 0;
            int b = primes[k];
            for (int m = i + 1; m > 0; m /= b) { f /= b; r += f * (m % b); }
            P[i][k] = std::fmod(r + shift[k], 1.0);
        }
    return P;
}

// ------------------------------------------------------------------ Gaussian process
double GaussianProcess::kernel(const std::vector<double>& a, const std::vector<double>& b) const {
    double r2 = 0;
    for (size_t k = 0; k < a.size(); k++) { double d = (a[k] - b[k]) / ls_[k]; r2 += d * d; }
    double r = std::sqrt(5.0 * r2);
    return sf2_ * (1 + r + r * r / 3) * std::exp(-r);
}

bool GaussianProcess::factor(const std::vector<double>& yn) {
    size_t n = X_.size();
    L_.assign(n * n, 0.0);
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j <= i; j++) L_[i * n + j] = kernel(X_[i], X_[j]) + (i == j ? sn2_ + 1e-10 * sf2_ : 0.0);
    for (size_t j = 0; j < n; j++) {
        double d = L_[j * n + j];
        for (size_t k = 0; k < j; k++) d -= L_[j * n + k] * L_[j * n + k];
        if (d <= 0) return false;
        d = std::sqrt(d);
        L_[j * n + j] = d;
        for (size_t i = j + 1; i < n; i++) {
            double s = L_[i * n + j];
            for (size_t k = 0; k < j; k++) s -= L_[i * n + k] * L_[j * n + k];
            L_[i * n + j] = s / d;
        }
    }
    alpha_ = yn;
    for (size_t i = 0; i < n; i++) { double s = alpha_[i]; for (size_t k = 0; k < i; k++) s -= L_[i * n + k] * alpha_[k]; alpha_[i] = s / L_[i * n + i]; }
    double fit = 0, logdet = 0;
    for (size_t i = 0; i < n; i++) { fit += alpha_[i] * alpha_[i]; logdet += std::log(L_[i * n + i]); }
    for (size_t i = n; i-- > 0;) { double s = alpha_[i]; for (size_t k = i + 1; k < n; k++) s -= L_[k * n + i] * alpha_[k]; alpha_[i] = s / L_[i * n + i]; }
    loglik_ = -0.5 * fit - logdet - 0.5 * n * std::log(2 * 3.14159265358979323846);
    return true;
}

static std::vector<double> normalised(const std::vector<double>& y, double& m, double& s) {
    m = std::accumulate(y.begin(), y.end(), 0.0) / y.size();
    double v = 0;
    for (double t : y) v += (t - m) * (t - m);
    s = std::sqrt(v / y.size());
    if (!(s > 0)) s = 1;
    std::vector<double> o(y.size());
    for (size_t i = 0; i < y.size(); i++) o[i] = (y[i] - m) / s;
    return o;
}

void GaussianProcess::fit(const std::vector<std::vector<double>>& X, const std::vector<double>& y, uint64_t seed) {
    if (X.size() != y.size() || X.empty()) throw std::invalid_argument("GP needs data");
    X_ = X;
    const size_t d = X[0].size();
    std::vector<double> yn = normalised(y, ymean_, ystd_);
    // theta = (log sf2, log l_1..l_d, log sn2), bounds as in the Python version
    const double lo_sf = std::log(1e-3), hi_sf = std::log(1e3), lo_l = std::log(1e-2), hi_l = std::log(1e2), lo_n = std::log(1e-8), hi_n = std::log(1.0);   // noise up to the variance of the data: DEM results scatter
    auto clampv = [&](std::vector<double>& t) {
        t[0] = std::min(std::max(t[0], lo_sf), hi_sf);
        for (size_t k = 0; k < d; k++) t[1 + k] = std::min(std::max(t[1 + k], lo_l), hi_l);
        t[d + 1] = std::min(std::max(t[d + 1], lo_n), hi_n);
    };
    auto objective = [&](std::vector<double> t) {
        clampv(t);
        sf2_ = std::exp(t[0]); ls_.resize(d);
        for (size_t k = 0; k < d; k++) ls_[k] = std::exp(t[1 + k]);
        sn2_ = std::exp(t[d + 1]);
        return factor(yn) ? -loglik_ : 1e300;
    };
    Rng rng(seed + 3);
    std::vector<double> best;
    double fbest = 1e300;
    const size_t np = d + 2;
    for (int start = 0; start < 6; start++) {
        std::vector<double> t0(np);
        if (start == 0) { t0[0] = 0; for (size_t k = 0; k < d; k++) t0[1 + k] = std::log(0.3); t0[d + 1] = std::log(1e-3); }
        else { t0[0] = rng.uniform(lo_sf, hi_sf) * 0.3; for (size_t k = 0; k < d; k++) t0[1 + k] = rng.uniform(std::log(0.05), std::log(3.0)); t0[d + 1] = rng.uniform(std::log(1e-6), std::log(1e-2)); }
        // Nelder-Mead
        std::vector<std::vector<double>> S(np + 1, t0);
        for (size_t k = 0; k < np; k++) S[k + 1][k] += 0.7;
        std::vector<double> fs(np + 1);
        for (size_t k = 0; k <= np; k++) fs[k] = objective(S[k]);
        for (int it = 0; it < 200 * int(np); it++) {
            std::vector<size_t> o(np + 1);
            std::iota(o.begin(), o.end(), 0);
            std::sort(o.begin(), o.end(), [&](size_t a, size_t b) { return fs[a] < fs[b]; });
            std::vector<std::vector<double>> S2; std::vector<double> f2;
            for (auto k : o) { S2.push_back(S[k]); f2.push_back(fs[k]); }
            S.swap(S2); fs.swap(f2);
            if (std::fabs(fs[np] - fs[0]) < 1e-7 * (1 + std::fabs(fs[0]))) break;
            std::vector<double> c(np, 0.0);
            for (size_t k = 0; k < np; k++) for (size_t j = 0; j < np; j++) c[j] += S[k][j] / np;
            auto pt = [&](double a) { std::vector<double> p(np); for (size_t j = 0; j < np; j++) p[j] = c[j] + a * (S[np][j] - c[j]); return p; };
            auto xr = pt(-1); double fr = objective(xr);
            if (fr < fs[0]) { auto xe = pt(-2); double fe = objective(xe); if (fe < fr) { S[np] = xe; fs[np] = fe; } else { S[np] = xr; fs[np] = fr; } }
            else if (fr < fs[np - 1]) { S[np] = xr; fs[np] = fr; }
            else {
                auto xc = pt(fr < fs[np] ? -0.5 : 0.5); double fc = objective(xc);
                if (fc < std::min(fr, fs[np])) { S[np] = xc; fs[np] = fc; }
                else for (size_t k = 1; k <= np; k++) { for (size_t j = 0; j < np; j++) S[k][j] = S[0][j] + 0.5 * (S[k][j] - S[0][j]); fs[k] = objective(S[k]); }
            }
        }
        if (fs[0] < fbest) { fbest = fs[0]; best = S[0]; }
    }
    objective(best);
}

void GaussianProcess::refit(const std::vector<std::vector<double>>& X, const std::vector<double>& y) {
    X_ = X;
    std::vector<double> yn(y.size());
    for (size_t i = 0; i < y.size(); i++) yn[i] = (y[i] - ymean_) / ystd_;
    if (!factor(yn)) throw std::runtime_error("GP refit failed");
}

void GaussianProcess::predict(const std::vector<double>& x, double& mean, double& sd) const {
    size_t n = X_.size();
    std::vector<double> k(n);
    for (size_t i = 0; i < n; i++) k[i] = kernel(x, X_[i]);
    double m = 0;
    for (size_t i = 0; i < n; i++) m += k[i] * alpha_[i];
    std::vector<double> v = k;
    for (size_t i = 0; i < n; i++) { double s = v[i]; for (size_t j = 0; j < i; j++) s -= L_[i * n + j] * v[j]; v[i] = s / L_[i * n + i]; }
    double var = sf2_;
    for (double t : v) var -= t * t;
    mean = ymean_ + ystd_ * m;
    sd = ystd_ * std::sqrt(std::max(var, 1e-16));
}

// ------------------------------------------------------------------ calibration
Calibration::Calibration(Evaluate f, std::vector<std::string> names, std::vector<double> lower, std::vector<double> upper)
    : f_(std::move(f)), names_(std::move(names)), lo_(std::move(lower)), hi_(std::move(upper)) {
    if (names_.empty() || names_.size() != lo_.size() || names_.size() != hi_.size()) throw std::invalid_argument("parameter names and bounds differ in size");
    for (size_t k = 0; k < lo_.size(); k++) if (!(hi_[k] > lo_[k])) throw std::invalid_argument("upper bound must exceed lower bound for " + names_[k]);
}
void Calibration::set_targets(const std::map<std::string, std::pair<double, double>>& t) {
    for (auto& kv : t) if (kv.second.first == 0) throw std::invalid_argument("target value of " + kv.first + " must not be zero");
    targets_ = t; objective_.clear();
    rescore();
}
void Calibration::set_objective(const std::string& o, const std::map<std::string, std::pair<double, double>>& c) { objective_ = o; constraints_ = c; targets_.clear(); rescore(); }

std::vector<double> Calibration::to_unit(const std::vector<double>& x) const {
    std::vector<double> u(x.size());
    for (size_t k = 0; k < x.size(); k++) u[k] = (x[k] - lo_[k]) / (hi_[k] - lo_[k]);
    return u;
}
std::vector<double> Calibration::from_unit(const std::vector<double>& u) const {
    std::vector<double> x(u.size());
    for (size_t k = 0; k < u.size(); k++) x[k] = lo_[k] + u[k] * (hi_[k] - lo_[k]);
    return x;
}

static double need(const std::map<std::string, double>& r, const std::string& k) {
    auto it = r.find(k);
    if (it == r.end()) throw std::runtime_error("quantity " + k + " missing in the evaluation result");
    return it->second;
}
double Calibration::misfit(const std::map<std::string, double>& r) const {
    double s = 0;
    for (auto& kv : targets_) { double e = (need(r, kv.first) - kv.second.first) / kv.second.first; s += kv.second.second * e * e; }
    return s;
}
bool Calibration::feasible(const std::map<std::string, double>& r) const {
    for (auto& kv : constraints_) { double v = need(r, kv.first); if (v < kv.second.first || v > kv.second.second) return false; }
    return true;
}
double Calibration::score(const std::map<std::string, double>& r) const {
    if (!targets_.empty()) return misfit(r);
    if (objective_.empty()) throw std::logic_error("set targets or an objective first");
    return need(r, objective_);
}

// Log line: label source x_1 ... x_d | key=value key=value ...
void Calibration::save(const Evaluation& e) const {
    if (log_.empty()) return;
    std::FILE* f = std::fopen(log_.c_str(), "a");
    if (!f) throw std::runtime_error("cannot write " + log_);
    std::fprintf(f, "%s %s", e.label.c_str(), e.source.c_str());
    for (double v : e.x) std::fprintf(f, " %.17g", v);
    std::fprintf(f, " |");
    for (auto& kv : e.result) std::fprintf(f, " %s=%.17g", kv.first.c_str(), kv.second);
    std::fprintf(f, "\n");
    std::fclose(f);
}

void Calibration::set_log(const std::string& path) {
    log_ = path;
    std::ifstream f(path);
    if (!f) return;
    std::string line;
    while (std::getline(f, line)) {
        auto bar = line.find('|');
        if (bar == std::string::npos) continue;
        std::istringstream a(line.substr(0, bar));
        Evaluation e;
        a >> e.label >> e.source;
        double v;
        while (a >> v) e.x.push_back(v);
        if (e.x.size() != names_.size()) throw std::runtime_error("log " + path + " has a different number of parameters");
        std::istringstream b(line.substr(bar + 1));
        std::string kv;
        while (b >> kv) { auto eq = kv.find('='); if (eq != std::string::npos) e.result[kv.substr(0, eq)] = std::stod(kv.substr(eq + 1)); }
        data_.push_back(e);
    }
    rescore();
}

void Calibration::rescore() {
    // scores depend on the targets or constraints, which may be set before or after the log is read
    if (targets_.empty() && objective_.empty()) return;
    for (auto& e : data_) {
        bool have = true;
        for (auto& kv : targets_) have &= e.result.count(kv.first) > 0;
        if (!objective_.empty()) have &= e.result.count(objective_) > 0;
        for (auto& kv : constraints_) have &= e.result.count(kv.first) > 0;
        if (!have) continue;
        e.score = score(e.result);
        e.feasible = feasible(e.result);
    }
}

void Calibration::evaluate(const std::vector<std::vector<double>>& X, const std::string& src) {
    for (auto& x : X) {
        Evaluation e;
        e.x = x;
        char lab[32];
        std::snprintf(lab, sizeof lab, "cal_%03zu", data_.size());
        e.label = lab; e.source = src;
        e.result = f_(x, e.label);
        e.score = score(e.result); e.feasible = feasible(e.result);
        data_.push_back(e);
        save(e);
    }
}

void Calibration::add(const std::vector<double>& x, const std::map<std::string, double>& r, const std::string& src) {
    Evaluation e;
    e.x = x; e.label = src; e.source = src; e.result = r; e.score = score(r); e.feasible = feasible(r);
    data_.push_back(e);
    save(e);
}

void Calibration::initial(int n) {
    for (auto& e : data_) { e.score = score(e.result); e.feasible = feasible(e.result); }
    if (!data_.empty()) return;
    auto U = latin_hypercube(n, int(names_.size()), seed_);
    std::vector<std::vector<double>> X;
    for (auto& u : U) X.push_back(from_unit(u));
    evaluate(X, "lhs");
}

std::vector<std::vector<double>> Calibration::propose(int it) {
    if (data_.size() < 2) throw std::logic_error("at least two evaluations are needed");
    const int d = int(names_.size());
    std::vector<std::vector<double>> U;
    for (auto& e : data_) U.push_back(to_unit(e.x));
    auto cand = halton(8192, d, seed_ + 100 + uint64_t(it));
    std::vector<std::vector<double>> picks;
    if (!targets_.empty()) {
        // One Gaussian process per measured quantity (with a noise term, so scattered DEM results are
        // averaged). For a candidate x the misfit S = sum w ((Y - y*)/y*)^2 with Y ~ N(m, s^2) has
        // mean sum w ((m - y*)^2 + s^2)/y*^2 and variance sum w^2 (2 s^4 + 4 s^2 (m - y*)^2)/y*^4;
        // S is approximated as normal and the expected improvement is taken over the smallest
        // expected misfit at the evaluated points (plug-in criterion for noisy observations).
        std::vector<std::string> q;
        std::vector<double> ts, ws;
        for (auto& kv : targets_) { q.push_back(kv.first); ts.push_back(kv.second.first); ws.push_back(kv.second.second); }
        std::vector<GaussianProcess> gp(q.size());
        std::vector<std::vector<double>> ys(q.size());
        for (size_t k = 0; k < q.size(); k++) {
            for (auto& e : data_) ys[k].push_back(need(e.result, q[k]));
            gp[k].fit(U, ys[k], seed_);
        }
        auto Ua = U;
        auto misfit = [&](const std::vector<double>& u, double& mean, double& sd) {
            mean = 0; double var = 0;
            for (size_t k = 0; k < q.size(); k++) {
                double m, s; gp[k].predict(u, m, s);
                double dlt = m - ts[k], t2 = ts[k] * ts[k];
                mean += ws[k] * (dlt * dlt + s * s) / t2;
                var += ws[k] * ws[k] * (2 * s * s * s * s + 4 * s * s * dlt * dlt) / (t2 * t2);
            }
            sd = std::sqrt(var);
        };
        for (int b = 0; b < batch_; b++) {
            double plug = std::numeric_limits<double>::infinity();
            for (auto& u : Ua) { double m, s; misfit(u, m, s); plug = std::min(plug, m); }
            double best = -1;
            size_t bi = 0;
            for (size_t i = 0; i < cand.size(); i++) {
                double m, s; misfit(cand[i], m, s);
                double a = expected_improvement(m, s, plug);
                if (a > best) { best = a; bi = i; }
            }
            picks.push_back(cand[bi]);
            // kriging believer: the predicted means are added as observations
            Ua.push_back(cand[bi]);
            for (size_t k = 0; k < q.size(); k++) {
                double m, s; gp[k].predict(cand[bi], m, s);
                ys[k].push_back(m);
                gp[k].refit(Ua, ys[k]);
            }
        }
    } else {
        std::vector<double> Y;
        bool anyfeas = false;
        for (auto& e : data_) { Y.push_back(need(e.result, objective_)); anyfeas |= e.feasible; }
        GaussianProcess go;
        go.fit(U, Y, seed_);
        std::vector<GaussianProcess> gc;
        std::vector<std::pair<double, double>> lims;
        for (auto& kv : constraints_) {
            std::vector<double> yk;
            for (auto& e : data_) yk.push_back(need(e.result, kv.first));
            gc.emplace_back(); gc.back().fit(U, yk, seed_);
            lims.push_back(kv.second);
        }
        double best = -std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < data_.size(); i++) if (data_[i].feasible) best = std::max(best, -Y[i]);
        best = anyfeas ? -best : *std::max_element(Y.begin(), Y.end());
        std::vector<double> acq(cand.size());
        for (size_t i = 0; i < cand.size(); i++) {
            double pf = 1, m, s;
            for (size_t c = 0; c < gc.size(); c++) {
                gc[c].predict(cand[i], m, s);
                s = std::max(s, 1e-12);
                pf *= normal_cdf((lims[c].second - m) / s) - normal_cdf((lims[c].first - m) / s);
            }
            go.predict(cand[i], m, s);
            acq[i] = anyfeas ? expected_improvement(m, s, best) * pf : pf;
        }
        for (int b = 0; b < batch_; b++) {
            size_t bi = size_t(std::max_element(acq.begin(), acq.end()) - acq.begin());
            picks.push_back(cand[bi]);
            for (size_t i = 0; i < cand.size(); i++) {
                double dd = 0;
                for (int k = 0; k < d; k++) dd += (cand[i][k] - cand[bi][k]) * (cand[i][k] - cand[bi][k]);
                if (std::sqrt(dd) < 0.15) acq[i] = -1;
            }
        }
    }
    std::vector<std::vector<double>> X;
    for (auto& u : picks) X.push_back(from_unit(u));
    return X;
}

void Calibration::step(int it) { evaluate(propose(it), "bo" + std::to_string(it)); }

Evaluation Calibration::run(int n_init, int n_iter, bool verbose) {
    initial(n_init);
    int start = 0;
    for (auto& e : data_) if (e.source.rfind("bo", 0) == 0) start++;
    start /= batch_;
    for (int it = start; it < n_iter; it++) {   // n_iter is the total: a resumed calibration completes it
        step(it);
        if (verbose) {
            Evaluation b = best();
            std::fprintf(stderr, "iteration %d: best", it);
            for (size_t k = 0; k < names_.size(); k++) std::fprintf(stderr, " %s=%.4g", names_[k].c_str(), b.x[k]);
            std::fprintf(stderr, "  score %.4g\n", b.score);
        }
    }
    return best();
}

Evaluation Calibration::best() const {
    if (data_.empty()) throw std::logic_error("no evaluations");
    const Evaluation* b = nullptr;
    bool any = false;
    for (auto& e : data_) any |= e.feasible;
    for (auto& e : data_) {
        if (!constraints_.empty() && any && !e.feasible) continue;
        if (!b || e.score < b->score) b = &e;
    }
    return *b;
}

std::map<std::string, std::vector<std::pair<double, double>>> Calibration::predict(const std::vector<std::vector<double>>& X) const {
    std::vector<std::vector<double>> U;
    for (auto& e : data_) U.push_back(to_unit(e.x));
    std::vector<std::string> q;
    for (auto& kv : targets_) q.push_back(kv.first);
    if (!objective_.empty()) q.push_back(objective_);
    for (auto& kv : constraints_) q.push_back(kv.first);
    std::map<std::string, std::vector<std::pair<double, double>>> out;
    for (auto& name : q) {
        std::vector<double> y;
        for (auto& e : data_) y.push_back(need(e.result, name));
        GaussianProcess gp;
        gp.fit(U, y, seed_);
        auto& v = out[name];
        for (auto& x : X) { double m, sd; gp.predict(to_unit(x), m, sd); v.push_back({m, sd}); }
    }
    return out;
}

Evaluation Calibration::estimate() const {
    // With scattered (noisy) results the best single run is often a lucky one. The estimate is the
    // minimum of the expected misfit (or of the objective) computed from the Gaussian processes of
    // the individual quantities, which average the scatter of neighbouring runs.
    if (data_.size() < 4) return best();
    const int d = int(names_.size());
    std::vector<std::vector<double>> U;
    for (auto& e : data_) U.push_back(to_unit(e.x));
    std::vector<std::string> q;
    if (!targets_.empty()) for (auto& kv : targets_) q.push_back(kv.first);
    else { q.push_back(objective_); for (auto& kv : constraints_) q.push_back(kv.first); }
    std::vector<GaussianProcess> gp(q.size());
    for (size_t k = 0; k < q.size(); k++) {
        std::vector<double> y;
        for (auto& e : data_) y.push_back(need(e.result, q[k]));
        gp[k].fit(U, y, seed_);
    }
    auto cand = halton(16384, d, seed_ + 7);
    for (auto& u : U) cand.push_back(u);
    double best_s = std::numeric_limits<double>::infinity();
    Evaluation out;
    for (auto& u : cand) {
        std::map<std::string, double> r;
        double unc = 0;   // expected misfit adds the surrogate variance: uncertain regions are not chosen
        for (size_t k = 0; k < q.size(); k++) {
            double m, sd; gp[k].predict(u, m, sd); r[q[k]] = m;
            auto t = targets_.find(q[k]);
            if (t != targets_.end()) unc += t->second.second * sd * sd / (t->second.first * t->second.first);
        }
        double sc;
        if (!targets_.empty()) sc = score(r) + unc;
        else { if (!feasible(r)) continue; sc = r[objective_]; }
        if (sc < best_s) { best_s = sc; out.x = from_unit(u); out.result = r; out.score = targets_.empty() ? score(r) : sc; }
    }
    if (out.x.empty()) return best();
    out.label = "estimate"; out.source = "gp";
    out.feasible = feasible(out.result);
    return out;
}

}  // namespace tb
