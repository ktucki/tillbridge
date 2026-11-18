// Minimal test framework (no external dependencies).
// Copyright (c) 2026 Karol Tucki. MIT License.
#pragma once
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

namespace tbt {
struct Case { const char* name; std::function<void()> fn; };
inline std::vector<Case>& cases() { static std::vector<Case> c; return c; }
inline int& failures() { static int f = 0; return f; }
struct Reg { Reg(const char* n, std::function<void()> f) { cases().push_back({n, std::move(f)}); } };
inline void fail(const char* file, int line, const std::string& msg) {
    std::printf("  FAILED %s:%d: %s\n", file, line, msg.c_str());
    failures()++;
}
}  // namespace tbt

#define TB_CAT2(a, b) a##b
#define TB_CAT(a, b) TB_CAT2(a, b)
#define TEST(name) \
    static void name(); \
    static tbt::Reg TB_CAT(reg_, name)(#name, name); \
    static void name()
#define CHECK(c) do { if (!(c)) tbt::fail(__FILE__, __LINE__, #c); } while (0)
#define CHECK_NEAR(a, b, tol) do { double _a = (a), _b = (b), _t = (tol); \
    if (!(std::fabs(_a - _b) <= _t)) tbt::fail(__FILE__, __LINE__, std::string(#a " = ") + std::to_string(_a) + ", expected " + std::to_string(_b) + " +- " + std::to_string(_t)); } while (0)
#define CHECK_REL(a, b, rtol) CHECK_NEAR(a, b, std::fabs(b) * (rtol))
