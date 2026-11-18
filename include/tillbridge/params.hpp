// TillBridge - solver parameters.
// Copyright (c) 2026 Karol Tucki. MIT License.
#pragma once
#include <string>
#include <vector>

namespace tb {

// Parameters of a DEM run (SI units). Defaults as in docs/file_formats.md.
struct Params {
    // material
    double E = 2.6e6, nu = 0.3, mu_pp = 0.6, mu_pw = 0.5, e = 0.4, mu_r = 0.3, k_coh = 0.0, E_stone = -1;
    double E_wall = 2.1e11, nu_wall = 0.3, g = 9.81;
    // run control
    double dt = 1e-4, t_end = 1.0, t_tool = 0.1, v_tool = 1.5, out_dt = 0.01;
    double Lx = 1.0, Ly = 1.0, Lz = 0.3, avg_x0 = 1e9, avg_x1 = -1e9, v_max = 50.0, ke_growth = 0.0, snap_dt = 0.0;
    // parallel bonds (Potyondy and Cundall 2004); bond_E = 0 switches them off
    double bond_E = 0, bond_kr = 2.5, bond_lam = 0.5, bond_sc = 1e30, bond_tc = 1e30, bond_gap = 0.05, bond_zeta = 0.2,
           bond_zmin = -1e9;
    std::string out = "tbdem";

    // Set one parameter by name; throws std::invalid_argument for an unknown key.
    void set(const std::string& key, const std::string& value);
    void set(const std::string& key, double value);
    double get(const std::string& key) const;
    // Read "key = value" lines ('#' starts a comment).
    static Params read(const std::string& path);
    void write(const std::string& path) const;
    static const std::vector<std::string>& keys();
    void check() const;
};

}  // namespace tb
