// Copyright (c) 2026 Karol Tucki. MIT License.
#include "tillbridge/params.hpp"

#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace tb {
namespace {

using Member = double Params::*;
const std::map<std::string, Member>& table() {
    static const std::map<std::string, Member> t = {
        {"E", &Params::E}, {"nu", &Params::nu}, {"mu_pp", &Params::mu_pp}, {"mu_pw", &Params::mu_pw},
        {"e", &Params::e}, {"mu_r", &Params::mu_r}, {"k_coh", &Params::k_coh}, {"E_stone", &Params::E_stone},
        {"E_wall", &Params::E_wall}, {"nu_wall", &Params::nu_wall}, {"g", &Params::g}, {"dt", &Params::dt},
        {"t_end", &Params::t_end}, {"t_tool", &Params::t_tool}, {"v_tool", &Params::v_tool},
        {"out_dt", &Params::out_dt}, {"Lx", &Params::Lx}, {"Ly", &Params::Ly}, {"Lz", &Params::Lz},
        {"avg_x0", &Params::avg_x0}, {"avg_x1", &Params::avg_x1}, {"v_max", &Params::v_max},
        {"ke_growth", &Params::ke_growth}, {"snap_dt", &Params::snap_dt}, {"bond_E", &Params::bond_E},
        {"bond_kr", &Params::bond_kr}, {"bond_lam", &Params::bond_lam}, {"bond_sc", &Params::bond_sc},
        {"bond_tc", &Params::bond_tc}, {"bond_gap", &Params::bond_gap}, {"bond_zeta", &Params::bond_zeta},
        {"bond_zmin", &Params::bond_zmin}};
    return t;
}

std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n"), b = s.find_last_not_of(" \t\r\n");
    return a == std::string::npos ? "" : s.substr(a, b - a + 1);
}

}  // namespace

const std::vector<std::string>& Params::keys() {
    static std::vector<std::string> k;
    if (k.empty())
        for (auto& p : table()) k.push_back(p.first);
    return k;
}

void Params::set(const std::string& key, double value) {
    auto it = table().find(key);
    if (it == table().end()) throw std::invalid_argument("unknown parameter " + key);
    this->*(it->second) = value;
}

void Params::set(const std::string& key, const std::string& value) {
    if (key == "out") { out = value; return; }
    size_t pos = 0;
    double d;
    try { d = std::stod(value, &pos); } catch (...) { throw std::invalid_argument("bad value for " + key + ": " + value); }
    if (trim(value.substr(pos)).size()) throw std::invalid_argument("bad value for " + key + ": " + value);
    set(key, d);
}

double Params::get(const std::string& key) const {
    auto it = table().find(key);
    if (it == table().end()) throw std::invalid_argument("unknown parameter " + key);
    return this->*(it->second);
}

Params Params::read(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open parameter file " + path);
    Params p;
    std::string line;
    while (std::getline(f, line)) {
        auto c = line.find('#');
        if (c != std::string::npos) line.resize(c);
        auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string k = trim(line.substr(0, eq)), v = trim(line.substr(eq + 1));
        if (k.empty()) continue;
        if (v.empty()) throw std::invalid_argument("missing value for " + k);
        // only the first token of the value counts (as in the C solver)
        std::istringstream is(v); std::string tok; is >> tok;
        p.set(k, tok);
    }
    p.check();
    return p;
}

void Params::write(const std::string& path) const {
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) throw std::runtime_error("cannot write " + path);
    for (auto& kv : table()) std::fprintf(f, "%s = %.17g\n", kv.first.c_str(), this->*(kv.second));
    std::fprintf(f, "out = %s\n", out.c_str());
    std::fclose(f);
}

void Params::check() const {
    if (dt <= 0 || t_end <= 0) throw std::invalid_argument("dt and t_end must be positive");
    if (E <= 0 || E_wall <= 0) throw std::invalid_argument("Young's moduli must be positive");
    if (nu < 0 || nu >= 0.5 || nu_wall < 0 || nu_wall >= 0.5) throw std::invalid_argument("Poisson's ratio out of range");
    if (Lx <= 0 || Ly <= 0) throw std::invalid_argument("box size must be positive");
}

}  // namespace tb
