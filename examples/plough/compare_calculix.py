"""Comparison of the TillBridge solution of the shank with CalculiX on the same mesh and nodal loads.

usage: python3 compare_calculix.py JOB VTK   (JOB.frd and JOB.dat from CalculiX, VTK from tillbridge fem --vtk)
Displacements are compared at the nodes; stresses as the largest von Mises stress at the Gauss points of
every element (same integration rules: 2x2x2 for C3D20R, 4 points for C3D10), and as nodal values
(both extrapolated from the Gauss points and averaged).
Copyright (c) 2026 Karol Tucki. MIT License.
"""
import sys, json
import numpy as np


def read_frd(fn):
    """Nodal displacements and stresses of the last step of an ASCII .frd file."""
    U, S = {}, {}
    block = None
    for line in open(fn, errors="replace"):
        if line.startswith(" -4"):
            name = line[5:13].strip()
            block = U if name == "DISP" else (S if name == "STRESS" else None)
            if block is not None: block.clear()
            continue
        if line.startswith(" -3"): block = None; continue
        if block is not None and line.startswith(" -1"):
            nid = int(line[3:13]); vals = [float(line[13 + 12 * k:25 + 12 * k]) for k in range((len(line.rstrip("\n")) - 13) // 12)]
            block[nid] = vals
    return U, S


def read_vtk(fn):
    L = open(fn).read().split("\n"); i = 0
    while not L[i].startswith("POINTS"): i += 1
    n = int(L[i].split()[1]); X = np.array([list(map(float, l.split())) for l in L[i + 1:i + 1 + n]])
    def field(key, ncomp):
        j = next(k for k, l in enumerate(L) if l.startswith(key))
        j += 2 if key.startswith("SCALARS") else 1
        return np.array([list(map(float, l.split())) for l in L[j:j + n]]).reshape(n, ncomp)
    j = next(k for k, l in enumerate(L) if l.startswith("CELL_DATA")); ne = int(L[j].split()[1])
    ve = np.array([float(l) for l in L[j + 3:j + 3 + ne]])
    return X, field("VECTORS", 3), field("SCALARS von_mises", 1)[:, 0], ve


def read_dat_stress(fn):
    """Largest von Mises stress over the integration points of every element (*EL PRINT S)."""
    out = {}; on = False
    for line in open(fn, errors="replace"):
        if "stresses (elem, integ.pnt." in line: on = True; continue
        if on:
            p = line.split()
            if len(p) == 8:
                e = int(p[0]); s = np.array(list(map(float, p[2:])))
                v = np.sqrt(0.5 * ((s[0] - s[1]) ** 2 + (s[1] - s[2]) ** 2 + (s[2] - s[0]) ** 2) + 3 * (s[3] ** 2 + s[4] ** 2 + s[5] ** 2))
                out[e] = max(out.get(e, 0.0), v)
            elif p: on = False
    return out


def vm(s):
    sx, sy, sz, txy, tyz, tzx = s.T
    return np.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2) + 3 * (txy ** 2 + tyz ** 2 + tzx ** 2))


case, vtk = sys.argv[1], sys.argv[2]
U, S = read_frd(f"{case}.frd")
X, Ut, VMt, VEt = read_vtk(vtk)
ED = read_dat_stress(f"{case}.dat")
VEc = np.array([ED[i] for i in range(1, len(VEt) + 1)])
ids = np.arange(1, len(X) + 1)
Uc = np.array([U[i][:3] for i in ids]); VMc = vm(np.array([S[i][:6] for i in ids]))
du = np.linalg.norm(Ut - Uc, axis=1)
tube = (X[:, 1] > 5) & (X[:, 1] < 420)
if not tube.any(): tube[:] = True
res = {
    "max_displacement_tillbridge": float(np.linalg.norm(Ut, axis=1).max()),
    "max_displacement_calculix": float(np.linalg.norm(Uc, axis=1).max()),
    "displacement_rel_l2": float(np.linalg.norm(du) / np.linalg.norm(Uc)),
    "max_vm_tube_tillbridge": float(VMt[tube].max()), "max_vm_tube_calculix": float(VMc[tube].max()),
    "max_vm_tillbridge": float(VMt.max()), "max_vm_calculix": float(VMc.max()),
    "vm_rel_l2": float(np.linalg.norm(VMt - VMc) / np.linalg.norm(VMc)),
    "vm_tube_rel_l2": float(np.linalg.norm(VMt[tube] - VMc[tube]) / np.linalg.norm(VMc[tube])),
    "element_vm_max_tillbridge": float(VEt.max()), "element_vm_max_calculix": float(VEc.max()),
    "element_vm_rel_l2": float(np.linalg.norm(VEt - VEc) / np.linalg.norm(VEc)),
    "element_vm_max_abs_diff": float(np.abs(VEt - VEc).max()),
    "reaction_check": "see JOB.dat",
}
print(json.dumps(res, indent=1))
json.dump(res, open(f"compare_{case}.json", "w"), indent=1)
np.savetxt(f"compare_{case}_elements.csv", np.c_[np.arange(1, len(VEt) + 1), VEt, VEc], delimiter=",",
           header="element,vm max Gauss TillBridge,vm max Gauss CalculiX", comments="", fmt="%.6g")
np.savetxt(f"compare_{case}_nodes.csv", np.c_[ids, X, np.linalg.norm(Ut, axis=1), np.linalg.norm(Uc, axis=1), VMt, VMc],
           delimiter=",", header="node,x,y,z,|u| TillBridge,|u| CalculiX,vm TillBridge,vm CalculiX", comments="", fmt="%.6g")
