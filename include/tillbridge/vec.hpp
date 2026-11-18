// TillBridge - small fixed-size vector and quaternion types.
// Copyright (c) 2026 Karol Tucki. MIT License.
#pragma once
#include <cmath>

namespace tb {

struct Vec3 {
    double x = 0, y = 0, z = 0;
    Vec3() = default;
    constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    double& operator[](int i) { return i == 0 ? x : (i == 1 ? y : z); }
    double operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }
    Vec3& operator+=(const Vec3& b) { x += b.x; y += b.y; z += b.z; return *this; }
    Vec3& operator-=(const Vec3& b) { x -= b.x; y -= b.y; z -= b.z; return *this; }
    Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
};

inline Vec3 operator+(const Vec3& a, const Vec3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator-(const Vec3& a, const Vec3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 operator-(const Vec3& a) { return {-a.x, -a.y, -a.z}; }
inline Vec3 operator*(const Vec3& a, double s) { return {a.x * s, a.y * s, a.z * s}; }
inline Vec3 operator*(double s, const Vec3& a) { return {a.x * s, a.y * s, a.z * s}; }
inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

// Unit quaternion (w, x, y, z) for body orientation.
struct Quat {
    double w = 1, x = 0, y = 0, z = 0;
};

// Rotation of p by q: p + 2w (u x p) + u x (2 u x p).
inline Vec3 rotate(const Quat& q, const Vec3& p) {
    Vec3 u{q.x, q.y, q.z};
    Vec3 t = cross(u, p) * 2.0;
    return p + t * q.w + cross(u, t);
}

inline Quat operator*(const Quat& a, const Quat& b) {
    return {a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z, a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
            a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x, a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w};
}

inline void normalize(Quat& q) {
    double n = std::sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
    q.w /= n; q.x /= n; q.y /= n; q.z /= n;
}

constexpr double kPi = 3.14159265358979323846;

}  // namespace tb
