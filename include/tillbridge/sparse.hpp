// TillBridge - compressed sparse row matrices, dense Cholesky and preconditioned CG.
// Copyright (c) 2026 Karol Tucki. MIT License.
#pragma once
#include <cstdint>
#include <functional>
#include <vector>

namespace tb {

struct Triplet {
    int64_t i, j;
    double v;
};

struct CsrMatrix {
    int64_t rows = 0, cols = 0;
    std::vector<int64_t> ptr;   // rows + 1
    std::vector<int64_t> idx;
    std::vector<double> val;

    // Duplicates are summed; column indices sorted within rows.
    static CsrMatrix from_triplets(int64_t rows, int64_t cols, std::vector<Triplet>& t);
    void multiply(const std::vector<double>& x, std::vector<double>& y) const;   // y = A x
    std::vector<double> diagonal() const;
    CsrMatrix transpose() const;
    // Submatrix A[keep, keep] (keep: sorted indices).
    CsrMatrix principal(const std::vector<int64_t>& keep) const;
    size_t nnz() const { return val.size(); }
};

CsrMatrix multiply(const CsrMatrix& A, const CsrMatrix& B);

// Dense symmetric positive definite matrix factorised in place (lower Cholesky).
struct DenseCholesky {
    int64_t n = 0;
    std::vector<double> L;
    void factor(std::vector<double> A, int64_t n);   // A row-major n x n
    void solve(std::vector<double>& b) const;          // in place
};

struct CgResult {
    int iterations = 0;
    double residual = 0;   // relative
    bool converged = false;
};

// Preconditioned conjugate gradients; precond(r, z) computes z = M^-1 r.
CgResult pcg(const CsrMatrix& A, const std::vector<double>& b, std::vector<double>& x,
             const std::function<void(const std::vector<double>&, std::vector<double>&)>& precond, double rtol = 1e-8,
             int maxiter = 5000);

double dot(const std::vector<double>& a, const std::vector<double>& b);
double norm2(const std::vector<double>& a);

}  // namespace tb
