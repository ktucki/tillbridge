// TillBridge - triangulated rigid tools.
// Copyright (c) 2026 Karol Tucki. MIT License.
#pragma once
#include <array>
#include <string>
#include <vector>

#include "tillbridge/vec.hpp"

namespace tb {

using Triangle = std::array<Vec3, 3>;

// Rigid tool made of triangles (metres, tool frame). The tool frame moves along +x at the
// travel speed; at t = t_tool its origin coincides with the origin of the soil box.
// File format: NT, then NT lines "ax ay az bx by bz cx cy cz".
struct Tool {
    std::vector<Triangle> tri;

    size_t size() const { return tri.size(); }
    Vec3 normal(size_t i) const;
    double area(size_t i) const;
    Vec3 centroid(size_t i) const;
    void bounds(Vec3& lo, Vec3& hi) const;

    Tool scaled(double s) const;
    Tool translated(const Vec3& d) const;
    // p' = R p + t (R row-major).
    Tool transformed(const double R[3][3], const Vec3& t) const;
    // Split triangles until no edge is longer than max_edge.
    Tool refined(double max_edge) const;
    Tool& append(const Tool& o);

    static Tool read(const std::string& path);
    // ASCII or binary STL; unit converts to metres (0.001 for mm).
    static Tool read_stl(const std::string& path, double unit = 1.0);
    void write(const std::string& path) const;
    void write_stl(const std::string& path, const std::string& name = "tool") const;
};

// Flat inclined tine with side faces; tip at (x_tip, y_centre, 0.3 - depth).
Tool chisel_tine(double width = 0.05, double depth = 0.15, double rake_deg = 60.0, double height = 0.35,
                 double y_centre = 0.0, double x_tip = 0.0, int n = 8);

// Generic mouldboard plough body (share, cylindroid mouldboard and landside), throwing the
// slice to -y. A parametric teaching geometry, not a particular commercial body.
struct MouldboardOptions {
    double width = 0.35, depth = 0.20, share_angle_deg = 40.0, top_angle_deg = 72.0, height_factor = 1.6;
    double rise = 0.55, curl = 1.2, spread = 0.35, landside_length = 0.75, y_land = 0.0, x_point = 0.0;
    double surface_z = 0.30;
    int nu = 12, nv = 12;
};
Tool mouldboard_body(const MouldboardOptions& o = {});

}  // namespace tb
