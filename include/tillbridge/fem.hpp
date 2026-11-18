// TillBridge - linear elastic solid finite elements for the structural member.
// Copyright (c) 2026 Karol Tucki. MIT License.
//
// Elements: 8-node and 20-node hexahedra (full 3x3x3 or reduced 2x2x2 Gauss integration) and
// 10-node tetrahedra. The system is solved by conjugate gradients preconditioned with a
// two-level smoothed-aggregation algebraic multigrid that carries the six rigid-body modes of
// every aggregate (Vanek, Mandel and Brezina 1996).
#pragma once
#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "tillbridge/sparse.hpp"
#include "tillbridge/vec.hpp"

namespace tb {

enum class ElementKind { Hex8, Hex20, Hex20Reduced, Tet10 };
int nodes_per_element(ElementKind k);
const char* element_name(ElementKind k);

struct ElementGroup {
    ElementKind kind;
    std::vector<int64_t> conn;   // nodes_per_element * n_elements, zero-based
    size_t size() const { return conn.size() / nodes_per_element(kind); }
};

struct Mesh {
    std::vector<Vec3> X;
    std::vector<ElementGroup> groups;
    size_t n_nodes() const { return X.size(); }
    size_t n_elements() const;
    double volume() const;
};

// Structured box mesh with nx*ny*nz hexahedra (order 1: Hex8, 2: Hex20).
Mesh box_mesh(const Vec3& L, int nx, int ny, int nz, const Vec3& origin = {}, int order = 2);
// 20-node hexahedra from 8-node ones (mid-edge nodes shared).
Mesh to_hex20(const Mesh& m);

struct Solution {
    std::vector<Vec3> U;                         // nodal displacements
    std::vector<std::vector<double>> vm_elem;    // von Mises stress per element (max over Gauss points), per group
    std::vector<double> vm_node;                 // nodally averaged von Mises stress
    Vec3 reaction;                               // total reaction at the fixed nodes
    int iterations = 0;
    double residual = 0;
};

struct SolverOptions {
    double rtol = 1e-8;
    int maxiter = 5000;
    double aggregate_size = 0;   // 0: automatic
    bool amg = true;             // false: Jacobi preconditioner
    bool verbose = false;
};

// Static linear analysis; fixed: nodes with all three displacements zero; loads: nodal forces.
Solution solve_static(const Mesh& mesh, double E, double nu, const std::vector<int64_t>& fixed,
                      const std::map<int64_t, Vec3>& loads, const SolverOptions& o = {});

CsrMatrix assemble_stiffness(const Mesh& mesh, double E, double nu);

// Shape functions and their derivatives in natural coordinates (exposed for testing).
void shape_functions(ElementKind k, double xi, double eta, double zeta, double* N, double (*dN)[3]);

// VTK legacy file with displacements and von Mises stress (for ParaView).
void write_vtk(const std::string& path, const Mesh& mesh, const Solution& s);

// The same linear static model as a CalculiX input file (C3D8, C3D20R, C3D10; nodes numbered from 1;
// node sets FIXED and LOADED; nodal forces as *CLOAD), for an independent check of the solution.
void write_calculix(const std::string& path, const Mesh& mesh, double E, double nu, const std::vector<int64_t>& fixed,
                    const std::map<int64_t, Vec3>& loads);

}  // namespace tb
