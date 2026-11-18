// Copyright (c) 2026 Karol Tucki. MIT License.
#include "tillbridge/sparse.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace tb {

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0;
    const int64_t n = int64_t(a.size());
#pragma omp parallel for reduction(+ : s) if (n > 20000)
    for (int64_t i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}
double norm2(const std::vector<double>& a) { return std::sqrt(dot(a, a)); }

CsrMatrix CsrMatrix::from_triplets(int64_t rows, int64_t cols, std::vector<Triplet>& t) {
    std::sort(t.begin(), t.end(), [](const Triplet& a, const Triplet& b) { return a.i != b.i ? a.i < b.i : a.j < b.j; });
    CsrMatrix A;
    A.rows = rows; A.cols = cols;
    A.ptr.assign(rows + 1, 0);
    for (size_t k = 0; k < t.size();) {
        if (t[k].i < 0 || t[k].i >= rows || t[k].j < 0 || t[k].j >= cols) throw std::out_of_range("triplet index out of range");
        size_t l = k;
        double v = 0;
        while (l < t.size() && t[l].i == t[k].i && t[l].j == t[k].j) v += t[l++].v;
        A.idx.push_back(t[k].j); A.val.push_back(v); A.ptr[t[k].i + 1]++;
        k = l;
    }
    for (int64_t i = 0; i < rows; i++) A.ptr[i + 1] += A.ptr[i];
    return A;
}

void CsrMatrix::multiply(const std::vector<double>& x, std::vector<double>& y) const {
    y.assign(rows, 0.0);
#pragma omp parallel for schedule(static) if (rows > 5000)
    for (int64_t i = 0; i < rows; i++) {
        double s = 0;
        for (int64_t k = ptr[i]; k < ptr[i + 1]; k++) s += val[k] * x[idx[k]];
        y[i] = s;
    }
}

std::vector<double> CsrMatrix::diagonal() const {
    std::vector<double> d(rows, 0.0);
    for (int64_t i = 0; i < rows; i++)
        for (int64_t k = ptr[i]; k < ptr[i + 1]; k++) if (idx[k] == i) d[i] = val[k];
    return d;
}

CsrMatrix CsrMatrix::transpose() const {
    CsrMatrix T;
    T.rows = cols; T.cols = rows;
    T.ptr.assign(cols + 1, 0);
    for (auto j : idx) T.ptr[j + 1]++;
    for (int64_t i = 0; i < cols; i++) T.ptr[i + 1] += T.ptr[i];
    T.idx.resize(idx.size()); T.val.resize(val.size());
    std::vector<int64_t> pos(T.ptr.begin(), T.ptr.end() - 1);
    for (int64_t i = 0; i < rows; i++)
        for (int64_t k = ptr[i]; k < ptr[i + 1]; k++) { int64_t p = pos[idx[k]]++; T.idx[p] = i; T.val[p] = val[k]; }
    return T;
}

CsrMatrix CsrMatrix::principal(const std::vector<int64_t>& keep) const {
    std::vector<int64_t> map(cols, -1);
    for (size_t k = 0; k < keep.size(); k++) map[keep[k]] = int64_t(k);
    CsrMatrix B;
    B.rows = B.cols = int64_t(keep.size());
    B.ptr.assign(B.rows + 1, 0);
    for (size_t r = 0; r < keep.size(); r++) {
        int64_t i = keep[r];
        for (int64_t k = ptr[i]; k < ptr[i + 1]; k++)
            if (map[idx[k]] >= 0) { B.idx.push_back(map[idx[k]]); B.val.push_back(val[k]); }
        B.ptr[r + 1] = int64_t(B.idx.size());
    }
    return B;
}

CsrMatrix multiply(const CsrMatrix& A, const CsrMatrix& B) {
    if (A.cols != B.rows) throw std::invalid_argument("matrix sizes do not match");
    CsrMatrix C;
    C.rows = A.rows; C.cols = B.cols;
    C.ptr.assign(A.rows + 1, 0);
    std::vector<double> acc(B.cols, 0.0);
    std::vector<int64_t> mark(B.cols, -1), cols;
    for (int64_t i = 0; i < A.rows; i++) {
        cols.clear();
        for (int64_t k = A.ptr[i]; k < A.ptr[i + 1]; k++) {
            int64_t a = A.idx[k];
            double v = A.val[k];
            for (int64_t l = B.ptr[a]; l < B.ptr[a + 1]; l++) {
                int64_t j = B.idx[l];
                if (mark[j] != i) { mark[j] = i; acc[j] = 0; cols.push_back(j); }
                acc[j] += v * B.val[l];
            }
        }
        std::sort(cols.begin(), cols.end());
        for (auto j : cols) { C.idx.push_back(j); C.val.push_back(acc[j]); }
        C.ptr[i + 1] = int64_t(C.idx.size());
    }
    return C;
}

void DenseCholesky::factor(std::vector<double> A, int64_t n_) {
    n = n_;
    for (int64_t j = 0; j < n; j++) {
        double d = A[j * n + j];
        for (int64_t k = 0; k < j; k++) d -= A[j * n + k] * A[j * n + k];
        if (d <= 0) throw std::runtime_error("matrix not positive definite (Cholesky)");
        d = std::sqrt(d);
        A[j * n + j] = d;
#pragma omp parallel for if (n - j > 256)
        for (int64_t i = j + 1; i < n; i++) {
            double s = A[i * n + j];
            for (int64_t k = 0; k < j; k++) s -= A[i * n + k] * A[j * n + k];
            A[i * n + j] = s / d;
        }
    }
    L.swap(A);
}

void DenseCholesky::solve(std::vector<double>& b) const {
    for (int64_t i = 0; i < n; i++) {
        double s = b[i];
        for (int64_t k = 0; k < i; k++) s -= L[i * n + k] * b[k];
        b[i] = s / L[i * n + i];
    }
    for (int64_t i = n - 1; i >= 0; i--) {
        double s = b[i];
        for (int64_t k = i + 1; k < n; k++) s -= L[k * n + i] * b[k];
        b[i] = s / L[i * n + i];
    }
}

CgResult pcg(const CsrMatrix& A, const std::vector<double>& b, std::vector<double>& x,
             const std::function<void(const std::vector<double>&, std::vector<double>&)>& M, double rtol, int maxiter) {
    const int64_t n = A.rows;
    if (int64_t(x.size()) != n) x.assign(n, 0.0);
    std::vector<double> r(n), z(n), p(n), q(n);
    A.multiply(x, q);
    for (int64_t i = 0; i < n; i++) r[i] = b[i] - q[i];
    double bn = norm2(b);
    CgResult res;
    if (bn == 0) { x.assign(n, 0.0); res.converged = true; return res; }
    M(r, z);
    p = z;
    double rz = dot(r, z);
    for (int it = 1; it <= maxiter; it++) {
        A.multiply(p, q);
        double alpha = rz / dot(p, q);
        for (int64_t i = 0; i < n; i++) { x[i] += alpha * p[i]; r[i] -= alpha * q[i]; }
        res.iterations = it;
        res.residual = norm2(r) / bn;
        if (res.residual < rtol) { res.converged = true; break; }
        M(r, z);
        double rz1 = dot(r, z), beta = rz1 / rz;
        rz = rz1;
        for (int64_t i = 0; i < n; i++) p[i] = z[i] + beta * p[i];
    }
    return res;
}

}  // namespace tb
