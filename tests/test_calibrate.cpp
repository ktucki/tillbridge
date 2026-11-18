// Gaussian-process calibration on analytic models with known answers.
// Copyright (c) 2026 Karol Tucki. MIT License.
#include <cmath>
#include <cstdio>

#include "testing.hpp"
#include "tillbridge/calibrate.hpp"

using namespace tb;

TEST(gp_interpolates_and_reports_uncertainty) {
    std::vector<std::vector<double>> X;
    std::vector<double> y;
    for (int i = 0; i <= 10; i++) { double x = i / 10.0; X.push_back({x}); y.push_back(std::sin(6 * x)); }
    GaussianProcess gp;
    gp.fit(X, y, 1);
    double m, s, m2, s2;
    gp.predict({0.55}, m, s);
    CHECK_NEAR(m, std::sin(3.3), 0.02);
    gp.predict({0.5}, m2, s2);
    CHECK(s2 < s + 1e-12);
    gp.predict({3.0}, m, s);
    CHECK(s > 10 * s2);
}

TEST(expected_improvement_properties) {
    CHECK_NEAR(normal_cdf(0), 0.5, 1e-15);
    CHECK_NEAR(normal_cdf(1.959963984540054), 0.975, 1e-12);
    CHECK_NEAR(expected_improvement(1.0, 0.0, 0.5), 0.0, 1e-15);
    CHECK_NEAR(expected_improvement(0.0, 0.0, 0.5), 0.5, 1e-15);
    CHECK_NEAR(expected_improvement(0.5, 1.0, 0.5), normal_pdf(0), 1e-15);
}

TEST(designs_fill_the_box) {
    auto L = latin_hypercube(10, 3, 7);
    for (int d = 0; d < 3; d++) {
        std::vector<int> bins(10, 0);
        for (auto& x : L) { CHECK(x[d] >= 0 && x[d] < 1); bins[int(x[d] * 10)]++; }
        for (int b : bins) CHECK(b == 1);
    }
    auto H = halton(64, 2, 3);
    double mean = 0;
    for (auto& x : H) mean += x[0] / 64;
    CHECK_NEAR(mean, 0.5, 0.03);
}

TEST(calibration_recovers_parameters) {
    // target matching recovers the parameters of a known smooth model
    auto model = [](const std::vector<double>& x, const std::string&) {
        return std::map<std::string, double>{{"draft", 4000 * x[0] + 0.010 * x[1]}, {"lift", 1500 + 2000 * x[0] - 0.006 * x[1]}};
    };
    std::vector<double> truth{0.55, 1.6e5};
    auto t = model(truth, "");
    Calibration cal(model, {"mu", "k"}, {0.2, 0.5e5}, {1.0, 2.5e5});
    cal.set_targets({{"draft", {t["draft"], 1.0}}, {"lift", {t["lift"], 1.0}}});
    cal.set_batch(2);
    cal.set_seed(1);
    Evaluation b = cal.run(8, 12, false);
    CHECK_NEAR(b.x[0], truth[0], 0.05);
    CHECK_NEAR(b.x[1] / truth[1], 1.0, 0.1);
}

TEST(constrained_calibration) {
    auto model = [](const std::vector<double>& x, const std::string&) {
        return std::map<std::string, double>{{"D", x[0] + x[1]}, {"C", std::pow(x[0] - 0.3, 2) + std::pow(x[1] - 0.7, 2)}};
    };
    Calibration cal(model, {"a", "b"}, {0, 0}, {1, 1});
    cal.set_objective("C", {{"D", {0.95, 1.05}}});
    cal.set_batch(2);
    Evaluation b = cal.run(8, 8, false);
    CHECK(b.feasible);
    CHECK(b.result["C"] < 0.01);
}

TEST(calibration_resumes_from_log) {
    int calls = 0;
    auto model = [&](const std::vector<double>& x, const std::string&) {
        calls++;
        return std::map<std::string, double>{{"y", x[0]}};
    };
    std::remove("tb_test_cal.log");
    {
        Calibration cal(model, {"a"}, {0}, {2});
        cal.set_targets({{"y", {1.2, 1.0}}});
        cal.set_log("tb_test_cal.log");
        cal.run(5, 3, false);
    }
    int first = calls;
    Calibration cal(model, {"a"}, {0}, {2});
    cal.set_log("tb_test_cal.log");          // log read before the targets are known
    cal.set_targets({{"y", {1.2, 1.0}}});
    CHECK(int(cal.data().size()) == first);
    cal.initial(5);   // skipped: data exist
    CHECK(calls == first);
    CHECK_NEAR(cal.best().x[0], 1.2, 0.1);
    CHECK(cal.best().score < 0.01);
    CHECK_NEAR(cal.estimate().x[0], 1.2, 0.05);
    std::remove("tb_test_cal.log");
}

TEST(estimate_averages_scattered_results) {
    // y = 100 (1 + 0.5 a) with +-5 % deterministic scatter; target 125 (a = 0.5). The best single run
    // can be a lucky one; the surrogate estimate must be close to the noise-free answer.
    auto model = [](const std::vector<double>& x, const std::string& label) {
        unsigned h = 2166136261u;
        for (char c : label) h = (h ^ unsigned(c)) * 16777619u;
        double noise = 0.05 * ((h % 2001) / 1000.0 - 1.0);
        return std::map<std::string, double>{{"y", 100 * (1 + 0.5 * x[0]) * (1 + noise)}};
    };
    Calibration cal(model, {"a"}, {0}, {1});
    cal.set_targets({{"y", {125.0, 1.0}}});
    cal.set_seed(5);
    cal.run(10, 10, false);
    CHECK_NEAR(cal.estimate().x[0], 0.5, 0.08);
}
