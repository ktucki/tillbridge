// Copyright (c) 2026 Karol Tucki. MIT License.
#include <chrono>
#include <cstring>

#include "testing.hpp"

int main(int argc, char** argv) {
    int run = 0;
    for (auto& c : tbt::cases()) {
        if (argc > 1 && !std::strstr(c.name, argv[1])) continue;
        int f0 = tbt::failures();
        auto t0 = std::chrono::steady_clock::now();
        try { c.fn(); } catch (std::exception& e) { tbt::fail(c.name, 0, std::string("exception: ") + e.what()); }
        double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::printf("%-44s %s  (%.2f s)\n", c.name, tbt::failures() == f0 ? "ok" : "FAILED", s);
        run++;
    }
    std::printf("\n%d tests, %d failed checks\n", run, tbt::failures());
    return tbt::failures() ? 1 : 0;
}
