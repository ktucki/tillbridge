// Copyright (c) 2026 Karol Tucki. MIT License.
#include "tillbridge/tool.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <fstream>
#include <functional>
#include <iterator>
#include <stdexcept>

#include "tillbridge/textio.hpp"

namespace tb {

Vec3 Tool::normal(size_t i) const {
    Vec3 n = cross(tri[i][1] - tri[i][0], tri[i][2] - tri[i][0]);
    double a = norm(n);
    return a > 0 ? n * (1.0 / a) : n;
}
double Tool::area(size_t i) const { return 0.5 * norm(cross(tri[i][1] - tri[i][0], tri[i][2] - tri[i][0])); }
Vec3 Tool::centroid(size_t i) const { return (tri[i][0] + tri[i][1] + tri[i][2]) * (1.0 / 3.0); }
void Tool::bounds(Vec3& lo, Vec3& hi) const {
    lo = {1e300, 1e300, 1e300}; hi = {-1e300, -1e300, -1e300};
    for (auto& t : tri)
        for (auto& p : t)
            for (int k = 0; k < 3; k++) { lo[k] = std::min(lo[k], p[k]); hi[k] = std::max(hi[k], p[k]); }
}

Tool Tool::scaled(double s) const { Tool o = *this; for (auto& t : o.tri) for (auto& p : t) p = p * s; return o; }
Tool Tool::translated(const Vec3& d) const { Tool o = *this; for (auto& t : o.tri) for (auto& p : t) p = p + d; return o; }
Tool Tool::transformed(const double R[3][3], const Vec3& t) const {
    Tool o = *this;
    for (auto& tr : o.tri)
        for (auto& p : tr) {
            Vec3 q{R[0][0] * p.x + R[0][1] * p.y + R[0][2] * p.z, R[1][0] * p.x + R[1][1] * p.y + R[1][2] * p.z,
                   R[2][0] * p.x + R[2][1] * p.y + R[2][2] * p.z};
            p = q + t;
        }
    return o;
}
Tool Tool::refined(double max_edge) const {
    if (max_edge <= 0) throw std::invalid_argument("max_edge must be positive");
    Tool o;
    std::vector<Triangle> stack(tri.rbegin(), tri.rend());
    while (!stack.empty()) {
        Triangle t = stack.back(); stack.pop_back();
        double e[3] = {norm(t[1] - t[0]), norm(t[2] - t[1]), norm(t[0] - t[2])};
        int k = int(std::max_element(e, e + 3) - e);
        if (e[k] <= max_edge) { o.tri.push_back(t); continue; }
        Vec3 a = t[k], b = t[(k + 1) % 3], c = t[(k + 2) % 3], m = (a + b) * 0.5;
        stack.push_back({m, b, c}); stack.push_back({a, m, c});
    }
    return o;
}
Tool& Tool::append(const Tool& o) { tri.insert(tri.end(), o.tri.begin(), o.tri.end()); return *this; }

Tool Tool::read(const std::string& path) {
    TextReader r(path);
    auto h = r.numbers();
    if (h.empty() || h[0] < 0) r.fail("bad tool header");
    Tool t; t.tri.resize(size_t(h[0]));
    for (auto& tr : t.tri) {
        auto v = r.numbers();
        if (v.size() < 9) r.fail("bad triangle record");
        for (int k = 0; k < 3; k++) tr[k] = {v[3 * k], v[3 * k + 1], v[3 * k + 2]};
    }
    return t;
}

Tool Tool::read_stl(const std::string& path, double unit) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::string data((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    Tool t;
    std::string head = data.substr(0, std::min<size_t>(400, data.size()));
    std::string h5 = data.substr(0, std::min<size_t>(5, data.size()));
    for (auto& c : h5) c = char(std::tolower(c));
    if (h5 == "solid" && head.find("facet") != std::string::npos) {
        std::vector<Vec3> v;
        size_t pos = 0;
        while ((pos = data.find("vertex", pos)) != std::string::npos) {
            pos += 6;
            auto nums = TextReader::parse(data.substr(pos, 120));
            if (nums.size() < 3) throw std::runtime_error("bad vertex in " + path);
            v.push_back({nums[0], nums[1], nums[2]});
        }
        if (v.size() % 3) throw std::runtime_error("bad facet count in " + path);
        for (size_t i = 0; i < v.size(); i += 3) t.tri.push_back({v[i] * unit, v[i + 1] * unit, v[i + 2] * unit});
        return t;
    }
    if (data.size() < 84) throw std::runtime_error("truncated binary STL " + path);
    uint32_t n; std::memcpy(&n, data.data() + 80, 4);
    if (data.size() < 84 + 50ull * n) throw std::runtime_error("truncated binary STL " + path);
    for (uint32_t i = 0; i < n; i++) {
        const char* p = data.data() + 84 + 50ull * i + 12;
        float v[9]; std::memcpy(v, p, 36);
        t.tri.push_back({Vec3{v[0], v[1], v[2]} * unit, Vec3{v[3], v[4], v[5]} * unit, Vec3{v[6], v[7], v[8]} * unit});
    }
    return t;
}

void Tool::write(const std::string& path) const {
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) throw std::runtime_error("cannot write " + path);
    std::fprintf(f, "%zu\n", tri.size());
    for (auto& t : tri)
        std::fprintf(f, "%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n", t[0].x, t[0].y, t[0].z, t[1].x, t[1].y, t[1].z,
                     t[2].x, t[2].y, t[2].z);
    std::fclose(f);
}

void Tool::write_stl(const std::string& path, const std::string& name) const {
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) throw std::runtime_error("cannot write " + path);
    std::fprintf(f, "solid %s\n", name.c_str());
    for (size_t i = 0; i < tri.size(); i++) {
        Vec3 n = normal(i);
        std::fprintf(f, " facet normal %.6e %.6e %.6e\n  outer loop\n", n.x, n.y, n.z);
        for (auto& p : tri[i]) std::fprintf(f, "   vertex %.6e %.6e %.6e\n", p.x, p.y, p.z);
        std::fprintf(f, "  endloop\n endfacet\n");
    }
    std::fprintf(f, "endsolid %s\n", name.c_str());
    std::fclose(f);
}

// Triangulate a parametric surface on an nu x nv grid of the unit square.
static Tool grid(const std::function<Vec3(double, double)>& fun, int nu, int nv) {
    Tool t;
    for (int i = 0; i < nu; i++)
        for (int j = 0; j < nv; j++) {
            double u0 = double(i) / nu, u1 = double(i + 1) / nu, v0 = double(j) / nv, v1 = double(j + 1) / nv;
            Vec3 a = fun(u0, v0), b = fun(u1, v0), c = fun(u1, v1), d = fun(u0, v1);
            t.tri.push_back({a, b, c}); t.tri.push_back({a, c, d});
        }
    return t;
}

Tool chisel_tine(double width, double depth, double rake_deg, double height, double y_centre, double x_tip, int n) {
    double th = rake_deg * kPi / 180.0, z0 = 0.3 - depth;
    Vec3 d{-std::cos(th), 0, std::sin(th)}, tip{x_tip, y_centre, z0};
    Tool t = grid([&](double u, double v) { return tip + d * (u * height) + Vec3{0, (v - 0.5) * width, 0}; }, n, 2);
    Vec3 thick = Vec3{-std::sin(th), 0, -std::cos(th)} * 0.02;
    for (double s : {-1.0, 1.0})
        t.append(grid([&](double u, double v) { return tip + d * (u * height) + Vec3{0, s * width / 2, 0} - thick * v; }, n, 1));
    return t;
}

Tool mouldboard_body(const MouldboardOptions& o) {
    double g0 = o.share_angle_deg * kPi / 180, g1 = o.top_angle_deg * kPi / 180;
    double H = o.height_factor * o.depth, zb = o.surface_z - o.depth;
    Tool t = grid([&](double U, double V) {
        double h = V * H, gam = g0 + (g1 - g0) * std::pow(V, 1.5);
        double xd = o.x_point - o.rise * h - o.curl * h * h;
        double s = U * o.width / std::sin(gam) * (1 + o.spread * V);
        return Vec3{xd - s * std::cos(gam), o.y_land - s * std::sin(gam), zb + h};
    }, o.nu, o.nv);
    double ztop = zb + 0.8 * H;
    t.append(grid([&](double U, double V) { return Vec3{o.x_point - U * o.landside_length, o.y_land, zb + V * (ztop - zb)}; }, 8, 3));
    return t;
}

}  // namespace tb
