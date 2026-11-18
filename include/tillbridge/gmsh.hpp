// TillBridge - meshes in the open Gmsh format (.msh, ASCII versions 2.2 and 4.1).
// Copyright (c) 2026 Karol Tucki. MIT License.
//
// Reading: 8-node and 20-node hexahedra, 10-node tetrahedra (Gmsh types 5, 17, 11), converted to the
// node order of the solver; every physical group, of any dimension, becomes a named node component
// (the nodes of its elements). Writing (version 2.2): solid elements and every node component as a
// physical group of one-node point elements, so a model read from ANSYS can be passed on to open tools.
#pragma once
#include <string>

#include "tillbridge/ansys.hpp"

namespace tb {

AnsysModel read_gmsh(const std::string& path);
void write_gmsh(const AnsysModel& model, const std::string& path);
// .msh -> Gmsh reader, anything else (.cdb, .dat) -> ANSYS reader.
AnsysModel read_model(const std::string& path);

}  // namespace tb
