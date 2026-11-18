// TillBridge - reading meshes written by ANSYS (CDBWRITE blocked format and Mechanical .dat).
// Copyright (c) 2026 Karol Tucki. MIT License.
//
// Only what a linear static model needs: nodes, solid elements (SOLID185/186/187 layouts with
// 8, 20 or 10 nodes), named node components, EX/NUXY constants and scalar parameters.
#pragma once
#include <map>
#include <string>
#include <vector>

#include "tillbridge/fem.hpp"

namespace tb {

struct AnsysModel {
    std::map<long, Vec3> nodes;
    std::vector<std::vector<long>> solids;                 // node ids of every solid element
    std::map<std::string, std::vector<long>> components;   // node components
    std::map<std::string, double> materials;               // e.g. EX, NUXY
    std::map<std::string, double> parameters;

    static AnsysModel read(const std::string& path);
    // Mesh with nodes renumbered 0..n-1; ids[k] is the ANSYS id of node k.
    Mesh mesh(std::vector<long>* ids = nullptr) const;
    // Zero-based node indices of a component.
    std::vector<int64_t> component(const std::string& name) const;
};

}  // namespace tb
