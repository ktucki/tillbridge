// Copyright (c) 2026 Karol Tucki. MIT License.
#include "tillbridge/gmsh.hpp"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace tb {
namespace {

// Node order of the solver (VTK / ANSYS) from the Gmsh order: solver[i] = gmsh[P[i]].
const int P_HEX20[20] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 9, 16, 18, 19, 17, 10, 12, 14, 15};
const int P_TET10[10] = {0, 1, 2, 3, 4, 5, 6, 7, 9, 8};

int solid_nodes(int type) {
    switch (type) {
        case 5: return 8;
        case 17: return 20;
        case 11: return 10;
        default: return 0;
    }
}
std::vector<long> to_solver(int type, const std::vector<long>& g) {
    std::vector<long> s(g.size());
    for (size_t i = 0; i < g.size(); i++) s[i] = type == 17 ? g[P_HEX20[i]] : (type == 11 ? g[P_TET10[i]] : g[i]);
    return s;
}
std::string lower(std::string s) { for (auto& c : s) c = char(std::tolower(static_cast<unsigned char>(c))); return s; }

struct Groups {
    std::map<std::pair<int, int>, std::string> names;            // (dim, physical tag) -> name
    std::map<std::pair<int, int>, std::set<long>> nodes;          // (dim, physical tag) -> node ids
    void add(int dim, int phys, const std::vector<long>& n) { nodes[{dim, phys}].insert(n.begin(), n.end()); }
    void store(AnsysModel& M) const {
        for (auto& g : nodes) {
            auto it = names.find(g.first);
            std::string name = it != names.end() ? it->second : "PHYS" + std::to_string(g.first.first) + "_" + std::to_string(g.first.second);
            auto& out = M.components[name];
            out.insert(out.end(), g.second.begin(), g.second.end());
            std::sort(out.begin(), out.end());
            out.erase(std::unique(out.begin(), out.end()), out.end());
        }
    }
};

void read_names(std::istream& f, Groups& G) {
    int n; f >> n;
    for (int k = 0; k < n; k++) {
        int dim, tag; f >> dim >> tag;
        std::string rest; std::getline(f, rest);
        size_t a = rest.find('"'), b = rest.rfind('"');
        G.names[{dim, tag}] = (a != std::string::npos && b > a) ? rest.substr(a + 1, b - a - 1) : std::to_string(tag);
    }
}

void read_v2(std::istream& f, AnsysModel& M, Groups& G) {
    std::string s;
    while (f >> s) {
        if (s == "$PhysicalNames") read_names(f, G);
        else if (s == "$Nodes") {
            long n; f >> n;
            for (long k = 0; k < n; k++) { long id; Vec3 x; f >> id >> x.x >> x.y >> x.z; M.nodes[id] = x; }
        } else if (s == "$Elements") {
            long n; f >> n;
            std::getline(f, s);
            for (long k = 0; k < n; k++) {
                std::getline(f, s);
                std::istringstream l(s);
                long id; int type, ntags; l >> id >> type >> ntags;
                std::vector<int> tags(ntags);
                for (auto& t : tags) l >> t;
                std::vector<long> nodes; for (long v; l >> v;) nodes.push_back(v);
                int dim = type == 15 ? 0 : (type == 1 || type == 8 ? 1 : (type == 2 || type == 3 || type == 9 || type == 16 || type == 10 ? 2 : 3));
                if (ntags > 0 && tags[0] != 0) G.add(dim, tags[0], nodes);
                int nn = solid_nodes(type);
                if (nn) {
                    if (int(nodes.size()) != nn) throw std::runtime_error("bad element " + std::to_string(id));
                    M.solids.push_back(to_solver(type, nodes));
                }
            }
        }
    }
}

void read_v4(std::istream& f, AnsysModel& M, Groups& G) {
    std::map<std::pair<int, int>, std::vector<int>> phys;   // (dim, entity) -> physical tags
    std::string s;
    while (f >> s) {
        if (s == "$PhysicalNames") read_names(f, G);
        else if (s == "$Entities") {
            long cnt[4]; f >> cnt[0] >> cnt[1] >> cnt[2] >> cnt[3];
            for (int dim = 0; dim < 4; dim++)
                for (long k = 0; k < cnt[dim]; k++) {
                    int tag; f >> tag;
                    double c; for (int j = 0; j < (dim == 0 ? 3 : 6); j++) f >> c;
                    int np; f >> np;
                    std::vector<int> p(np); for (auto& t : p) f >> t;
                    phys[{dim, tag}] = p;
                    if (dim > 0) { int nb; f >> nb; for (int j = 0; j < nb; j++) { int b; f >> b; } }
                }
        } else if (s == "$Nodes") {
            long blocks, n, mn, mx; f >> blocks >> n >> mn >> mx;
            for (long b = 0; b < blocks; b++) {
                int dim, ent, par; long nb; f >> dim >> ent >> par >> nb;
                std::vector<long> ids(nb); for (auto& t : ids) f >> t;
                for (long k = 0; k < nb; k++) {
                    Vec3 x; f >> x.x >> x.y >> x.z;
                    double u; for (int j = 0; j < (par ? dim : 0); j++) f >> u;
                    M.nodes[ids[k]] = x;
                }
            }
        } else if (s == "$Elements") {
            long blocks, n, mn, mx; f >> blocks >> n >> mn >> mx;
            std::getline(f, s);
            for (long b = 0; b < blocks; b++) {
                std::getline(f, s);
                std::istringstream h(s);
                int dim, ent, type; long nb; h >> dim >> ent >> type >> nb;
                const auto& p = phys[{dim, ent}];
                int nn = solid_nodes(type);
                for (long k = 0; k < nb; k++) {
                    std::getline(f, s);
                    std::istringstream l(s);
                    long id; l >> id;
                    std::vector<long> nodes; for (long v; l >> v;) nodes.push_back(v);
                    for (int t : p) G.add(dim, t, nodes);
                    if (nn) {
                        if (int(nodes.size()) != nn) throw std::runtime_error("bad element " + std::to_string(id));
                        M.solids.push_back(to_solver(type, nodes));
                    }
                }
            }
        }
    }
}

}  // namespace

AnsysModel read_gmsh(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::string s; double version = 0; int binary = 0, dsize;
    while (f >> s) if (s == "$MeshFormat") { f >> version >> binary >> dsize; break; }
    if (version == 0) throw std::runtime_error(path + " is not a Gmsh mesh");
    if (binary) throw std::runtime_error(path + ": binary .msh is not supported, save the mesh as ASCII");
    AnsysModel M; Groups G;
    if (version < 3) read_v2(f, M, G);
    else if (version >= 4) read_v4(f, M, G);
    else throw std::runtime_error(path + ": .msh version " + std::to_string(version) + " is not supported");
    G.store(M);
    if (M.solids.empty()) throw std::runtime_error(path + " has no hexahedra or tetrahedra of a supported type");
    return M;
}

void write_gmsh(const AnsysModel& M, const std::string& path) {
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) throw std::runtime_error("cannot write " + path);
    std::fprintf(f, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n%zu\n3 1 \"SOLID\"\n", M.components.size() + 1);
    int tag = 2;
    for (auto& c : M.components) std::fprintf(f, "0 %d \"%s\"\n", tag++, c.first.c_str());
    std::fprintf(f, "$EndPhysicalNames\n$Nodes\n%zu\n", M.nodes.size());
    for (auto& n : M.nodes) std::fprintf(f, "%ld %.10g %.10g %.10g\n", n.first, n.second.x, n.second.y, n.second.z);
    size_t ne = M.solids.size();
    for (auto& c : M.components) ne += c.second.size();
    std::fprintf(f, "$EndNodes\n$Elements\n%zu\n", ne);
    long id = 1;
    for (auto& e : M.solids) {
        int type = e.size() == 20 ? 17 : (e.size() == 10 ? 11 : 5);
        std::vector<long> g(e.size());
        for (size_t i = 0; i < e.size(); i++) {
            int j = type == 17 ? P_HEX20[i] : (type == 11 ? P_TET10[i] : int(i));
            g[j] = e[i];
        }
        std::fprintf(f, "%ld %d 2 1 1", id++, type);
        for (long v : g) std::fprintf(f, " %ld", v);
        std::fprintf(f, "\n");
    }
    tag = 2;
    for (auto& c : M.components) {
        for (long v : c.second) std::fprintf(f, "%ld 15 2 %d %d %ld\n", id++, tag, tag, v);
        tag++;
    }
    std::fprintf(f, "$EndElements\n");
    std::fclose(f);
}

AnsysModel read_model(const std::string& path) {
    std::string l = lower(path);
    if (l.size() > 4 && l.compare(l.size() - 4, 4, ".msh") == 0) return read_gmsh(path);
    return AnsysModel::read(path);
}

}  // namespace tb
