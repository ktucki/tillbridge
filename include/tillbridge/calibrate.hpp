// TillBridge - calibration of DEM parameters by Bayesian optimisation.
// Copyright (c) 2026 Karol Tucki. MIT License.
//
// Every simulation is expensive, so parameters are searched with Gaussian-process surrogates
// (Matern 5/2 kernel with one length scale per parameter, hyperparameters by maximum marginal
// likelihood). Two problem types:
//  * target matching: minimise sum w ((y - y*) / y*)^2 over measured quantities;
//  * constrained: minimise one quantity while others stay in intervals (expected improvement
//    times the probability of feasibility; Gardner et al. 2014, Gelbart et al. 2014).
// Batches are proposed with the kriging believer, so several runs can be evaluated together.
// Every evaluation is appended to a text log, from which a calibration can be resumed.
#pragma once
#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace tb {

// Gaussian process on the unit box.
class GaussianProcess {
  public:
    // Fits hyperparameters (log-likelihood maximised by Nelder-Mead from several starts).
    void fit(const std::vector<std::vector<double>>& X, const std::vector<double>& y, uint64_t seed = 0);
    // Refit with fixed hyperparameters (for the kriging believer).
    void refit(const std::vector<std::vector<double>>& X, const std::vector<double>& y);
    void predict(const std::vector<double>& x, double& mean, double& sd) const;
    double log_likelihood() const { return loglik_; }
    const std::vector<double>& length_scales() const { return ls_; }

  private:
    std::vector<std::vector<double>> X_;
    std::vector<double> alpha_, L_, ls_;
    double sf2_ = 1, sn2_ = 1e-3, ymean_ = 0, ystd_ = 1, loglik_ = 0;
    bool factor(const std::vector<double>& yn);
    double kernel(const std::vector<double>& a, const std::vector<double>& b) const;
};

double normal_pdf(double z);
double normal_cdf(double z);
double expected_improvement(double mean, double sd, double best);
// Latin hypercube design of n points in [0, 1]^d.
std::vector<std::vector<double>> latin_hypercube(int n, int d, uint64_t seed);
// Halton points with a random shift (Cranley-Patterson rotation).
std::vector<std::vector<double>> halton(int n, int d, uint64_t seed);

struct Evaluation {
    std::vector<double> x;
    std::string label, source;
    std::map<std::string, double> result;
    double score = 0;
    bool feasible = true;
};

class Calibration {
  public:
    using Evaluate = std::function<std::map<std::string, double>(const std::vector<double>& x, const std::string& label)>;

    Calibration(Evaluate f, std::vector<std::string> names, std::vector<double> lower, std::vector<double> upper);
    // Target matching: quantity -> (measured value, weight).
    void set_targets(const std::map<std::string, std::pair<double, double>>& t);
    // Constrained: minimise `objective` subject to quantity -> (low, high).
    void set_objective(const std::string& objective, const std::map<std::string, std::pair<double, double>>& constraints);
    void set_log(const std::string& path);   // resumes from an existing log
    void set_batch(int b) { batch_ = b < 1 ? 1 : b; }
    void set_seed(uint64_t s) { seed_ = s; }

    void initial(int n);                                    // Latin hypercube (skipped if data exist)
    std::vector<std::vector<double>> propose(int iteration);
    void step(int iteration);
    // n_init start points and n_iter Bayesian steps in total; a resumed calibration only adds
    // what is missing. Returns the best run; see estimate() for scattered results.
    Evaluation run(int n_init, int n_iter, bool verbose = true);
    Evaluation best() const;       // best evaluated run
    Evaluation estimate() const;   // minimum of the surrogate means (use with scattered results)
    // Surrogate mean and standard deviation of every target, objective or constraint quantity at
    // the points X (parameter units); one Gaussian process per quantity, fitted to all runs.
    std::map<std::string, std::vector<std::pair<double, double>>> predict(const std::vector<std::vector<double>>& X) const;
    void add(const std::vector<double>& x, const std::map<std::string, double>& result, const std::string& source = "user");
    const std::vector<Evaluation>& data() const { return data_; }
    double misfit(const std::map<std::string, double>& r) const;
    bool feasible(const std::map<std::string, double>& r) const;

  private:
    Evaluate f_;
    std::vector<std::string> names_;
    std::vector<double> lo_, hi_;
    std::map<std::string, std::pair<double, double>> targets_, constraints_;
    std::string objective_, log_;
    int batch_ = 1;
    uint64_t seed_ = 0;
    std::vector<Evaluation> data_;
    double score(const std::map<std::string, double>& r) const;
    void rescore();
    void evaluate(const std::vector<std::vector<double>>& X, const std::string& source);
    void save(const Evaluation& e) const;
    std::vector<double> to_unit(const std::vector<double>& x) const;
    std::vector<double> from_unit(const std::vector<double>& u) const;
};

}  // namespace tb
