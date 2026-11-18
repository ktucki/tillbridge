// Copyright (c) 2026 Karol Tucki. MIT License.
#include "tillbridge/ansys.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <regex>
#include <stdexcept>

namespace tb {
namespace {

std::string lower(std::string s) { for (auto& c : s) c = char(std::tolower(static_cast<unsigned char>(c))); return s; }
std::string upper(std::string s) { for (auto& c : s) c = char(std::toupper(static_cast<unsigned char>(c))); return s; }
std::string strip(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r"), b = s.find_last_not_of(" \t\r");
    return a == std::string::npos ? "" : s.substr(a, b - a + 1);
}

// Fixed-width integers of a record.
std::vector<long> ints(const std::string& r0, int w) {
    std::string r = r0;
    while (!r.empty() && (r.back() == ' ' || r.back() == '\r')) r.pop_back();
    std::vector<long> v;
    for (size_t k = 0; k < r.size(); k += w) {
        std::string f = strip(r.substr(k, w));
        if (!f.empty()) v.push_back(std::stol(f));
    }
    return v;
}

int width_i(const std::string& fmt, int def) {
    std::smatch m;
    std::string l = lower(fmt);
    if (std::regex_search(l, m, std::regex("i(\\d+)"))) return std::stoi(m[1]);
    return def;
}

}  // namespace

AnsysModel AnsysModel::read(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::vector<std::string> L;
    for (std::string s; std::getline(f, s);) { if (!s.empty() && s.back() == '\r') s.pop_back(); L.push_back(s); }
    AnsysModel M;
    const std::regex par("\\s*(\\w+)\\s*=\\s*([-+0-9.eE]+)\\s*$"), tab("(_\\w+)\\(2,1,1\\)\\s*=\\s*([-+0-9.eE]+).*");
    size_t i = 0;
    while (i < L.size()) {
        const std::string& s = L[i];
        std::string sl = lower(s);
        if (sl.rfind("nblock", 0) == 0) {
            std::string fmt = lower(L.at(i + 1));
            i += 2;
            int ni = 3, wi = 9, wf = 20;
            std::smatch m;
            if (std::regex_search(fmt, m, std::regex("(\\d+)i(\\d+)"))) { ni = std::stoi(m[1]); wi = std::stoi(m[2]); }
            if (std::regex_search(fmt, m, std::regex("(\\d+)e(\\d+)"))) wf = std::stoi(m[2]);
            while (i < L.size() && L[i].rfind("-1", 0) != 0 && upper(L[i]).rfind("N,", 0) != 0) {
                const std::string& r = L[i];
                if (!strip(r).empty()) {
                    long id = std::stol(r.substr(0, wi));
                    size_t off = size_t(ni) * wi;
                    Vec3 x;
                    for (int k = 0; k < 3; k++) {
                        std::string fld = off + k * wf < r.size() ? strip(r.substr(off + k * wf, wf)) : "";
                        x[k] = fld.empty() ? 0.0 : std::stod(fld);
                    }
                    M.nodes[id] = x;
                }
                i++;
            }
        } else if (sl.rfind("eblock", 0) == 0) {
            bool solid = sl.find("solid") != std::string::npos;
            int w = width_i(L.at(i + 1), 9);
            i += 2;
            while (i < L.size() && strip(L[i]).rfind("-1", 0) != 0) {
                auto v = ints(L[i], w);
                if (solid && v.size() >= 11) {
                    size_t nn = size_t(v[8]);
                    std::vector<long> conn(v.begin() + 11, v.end());
                    while (conn.size() < nn && i + 1 < L.size()) { i++; auto more = ints(L[i], w); conn.insert(conn.end(), more.begin(), more.end()); }
                    conn.resize(nn);
                    M.solids.push_back(conn);
                }
                i++;
            }
        } else if (sl.rfind("cmblock", 0) == 0 && sl.find("node") != std::string::npos) {
            std::vector<std::string> p;
            size_t a = 0, b;
            while ((b = s.find(',', a)) != std::string::npos) { p.push_back(s.substr(a, b - a)); a = b + 1; }
            p.push_back(s.substr(a));
            if (p.size() < 4) throw std::runtime_error("bad CMBLOCK in " + path);
            std::string name = strip(p[1]);
            size_t cnt = size_t(std::stol(strip(p[3])));
            int w = width_i(L.at(i + 1), 10);
            i += 2;
            std::vector<long> ids;
            while (ids.size() < cnt && i < L.size()) { auto v = ints(L[i], w); ids.insert(ids.end(), v.begin(), v.end()); i++; }
            std::vector<long> out;
            for (long x : ids) {
                if (x < 0 && !out.empty()) for (long k = out.back() + 1; k <= -x; k++) out.push_back(k);
                else out.push_back(x);
            }
            M.components[name] = out;
            continue;
        } else if (upper(s).rfind("MP,", 0) == 0) {
            std::vector<std::string> p;
            size_t a = 0, b;
            while ((b = s.find(',', a)) != std::string::npos) { p.push_back(s.substr(a, b - a)); a = b + 1; }
            p.push_back(s.substr(a));
            if (p.size() > 3) { try { M.materials[upper(strip(p[1]))] = std::stod(p[3]); } catch (...) {} }
        } else {
            std::smatch m;
            if (std::regex_match(s, m, par)) { try { M.parameters[m[1]] = std::stod(m[2]); } catch (...) {} }
            if (std::regex_match(s, m, tab)) { try { M.parameters[m[1]] = std::stod(m[2]); } catch (...) {} }
        }
        i++;
    }
    return M;
}

Mesh AnsysModel::mesh(std::vector<long>* ids_out) const {
    Mesh m;
    std::map<long, int64_t> idx;
    std::vector<long> ids;
    for (auto& n : nodes) { idx[n.first] = int64_t(m.X.size()); m.X.push_back(n.second); ids.push_back(n.first); }
    ElementGroup h20{ElementKind::Hex20Reduced, {}}, h8{ElementKind::Hex8, {}}, t10{ElementKind::Tet10, {}};
    for (auto& e : solids) {
        ElementGroup* g = e.size() == 20 ? &h20 : (e.size() == 8 ? &h8 : (e.size() == 10 ? &t10 : nullptr));
        if (!g) continue;
        for (long id : e) {
            auto it = idx.find(id);
            if (it == idx.end()) throw std::runtime_error("element refers to an unknown node " + std::to_string(id));
            g->conn.push_back(it->second);
        }
    }
    for (auto* g : {&h20, &h8, &t10}) if (!g->conn.empty()) m.groups.push_back(*g);
    if (ids_out) *ids_out = ids;
    return m;
}

std::vector<int64_t> AnsysModel::component(const std::string& name) const {
    auto it = components.find(name);
    if (it == components.end()) throw std::runtime_error("no node component " + name);
    std::map<long, int64_t> idx;
    int64_t k = 0;
    for (auto& n : nodes) idx[n.first] = k++;
    std::vector<int64_t> out;
    for (long id : it->second) { auto j = idx.find(id); if (j != idx.end()) out.push_back(j->second); }
    return out;
}

}  // namespace tb
