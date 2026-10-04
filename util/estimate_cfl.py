#!/usr/bin/env python3
"""
由 Mesh3d.dat (无格式) + flow3d.dat 估算 CFL 数
几何量、谱半径的计算与 sub_geometry.f90 / sub_time_acceleraction.f90 (comput_Lijk, comput_dt) 一致:
    CFL = dt * (Lci+Lcj+Lck + Lvi+Lvj+Lvk) / Vol
另外输出只含对流项的 CFL_c = dt * (Lci+Lcj+Lck) / Vol
注: 湍流粘性 mu_t 不在 flow3d.dat 中, 按 0 处理 (层流精确, 湍流算例粘性部分偏小)

用法:  python3 estimate_cfl.py [--dir 算例目录] [--dt 时间步长] [--flow flow3d.dat] [--freestream]
       dt 默认取 control.ec 中的 dt_global
"""
import argparse
import os
import re
import sys

import numpy as np


def read_namelist(fname):
    """简单解析 control.ec 中的 key=value"""
    para = {}
    with open(fname) as f:
        for line in f:
            line = line.split("!")[0]
            for key, val in re.findall(r"(\w+)\s*=\s*([^\s,/]+)", line):
                para[key.lower()] = val.lower().replace("d", "e")
    return para


class FortranFile:
    """顺序读取 Fortran 无格式文件 (4字节记录标记)"""

    def __init__(self, fname):
        self.f = open(fname, "rb")

    def read(self, dtype):
        n = np.fromfile(self.f, dtype="<i4", count=1)
        if n.size == 0:
            raise EOFError
        data = np.fromfile(self.f, dtype=dtype, count=n[0] // np.dtype(dtype).itemsize)
        n2 = np.fromfile(self.f, dtype="<i4", count=1)
        if n2.size == 0 or n2[0] != n[0]:
            raise IOError("记录标记不匹配, 文件格式错误或文件不完整")
        return data

    def close(self):
        self.f.close()


def block_geometry(x, y, z):
    """x,y,z: (nx,ny,nz) 节点坐标; 返回各方向面积、法向及单元体积"""
    P = np.stack([x, y, z], axis=-1)
    nx, ny, nz = x.shape
    I, J, K = slice(0, nx - 1), slice(0, ny - 1), slice(0, nz - 1)
    I1, J1, K1 = slice(1, nx), slice(1, ny), slice(1, nz)

    # i 面: i=1..nx, j=1..ny-1, k=1..nz-1  (与 sub_geometry.f90 中对角线叉乘方向一致)
    p00, p10 = P[:, J, K], P[:, J1, K]
    p01, p11 = P[:, J, K1], P[:, J1, K1]
    s1 = np.cross(p10 - p01, p11 - p00)
    ctr = 0.25 * (p00 + p10 + p01 + p11)
    Si, ni, Vi = _face(s1, ctr)

    # j 面: i=1..nx-1, j=1..ny, k=1..nz-1
    p00, p10 = P[I, :, K], P[I1, :, K]
    p01, p11 = P[I, :, K1], P[I1, :, K1]
    s1 = np.cross(p11 - p00, p10 - p01)
    ctr = 0.25 * (p00 + p10 + p01 + p11)
    Sj, nj, Vj = _face(s1, ctr)

    # k 面: i=1..nx-1, j=1..ny-1, k=1..nz
    p00, p10 = P[I, J, :], P[I1, J, :]
    p01, p11 = P[I, J1, :], P[I1, J1, :]
    s1 = np.cross(p11 - p00, p01 - p10)
    ctr = 0.25 * (p00 + p10 + p01 + p11)
    Sk, nk, Vk = _face(s1, ctr)

    vol = (Vi[1:, :, :] - Vi[:-1, :, :] + Vj[:, 1:, :] - Vj[:, :-1, :] + Vk[:, :, 1:] - Vk[:, :, :-1]) / 3.0
    return (Si, ni), (Sj, nj), (Sk, nk), vol


def _face(s1, ctr):
    ss = np.sqrt(np.sum(s1 * s1, axis=-1))
    with np.errstate(invalid="ignore", divide="ignore"):
        n = np.where(ss[..., None] >= 1.0e-30, s1 / ss[..., None], np.array([1.0, 0.0, 0.0]))
    V = 0.5 * np.sum(s1 * ctr, axis=-1)
    return 0.5 * ss, n, V


def spectral_radius(geo, vol, rho, u, v, w, T, para):
    """与 comput_Lijk 一致: 返回 (Lc, Lv) = (Lci+Lcj+Lck, Lvi+Lvj+Lvk)"""
    (Si, ni), (Sj, nj), (Sk, nk), _ = geo
    Ma, gamma, Re, PrL = para["ma"], para["gamma"], para["re"], para["prl"]
    vel = np.stack([u, v, w], axis=-1)
    cc = np.sqrt(T) / Ma
    if para["if_viscous"] == 1:
        Tsb = 110.4 / para["t_inf"]
        mu = 1.0 / Re * (1.0 + Tsb) * np.sqrt(T**3) / (Tsb + T)
        D0 = max(gamma, 4.0 / 3.0) / rho * (mu / PrL)  # mu_t 未知, 取 0
    else:
        D0 = np.zeros_like(rho)
    Lc = np.zeros_like(rho)
    Lv = np.zeros_like(rho)
    for S, n, ax in ((Si, ni, 0), (Sj, nj, 1), (Sk, nk, 2)):
        lo = [slice(None)] * 3
        hi = [slice(None)] * 3
        lo[ax] = slice(0, -1)
        hi[ax] = slice(1, None)
        nm = n[tuple(lo)] + n[tuple(hi)]
        un = 0.5 * np.sum(vel * nm, axis=-1)
        S0 = 0.5 * (S[tuple(lo)] + S[tuple(hi)])
        Lc += (np.abs(un) + cc) * S0
        Lv += D0 * S0 * S0 / vol
    return Lc, Lv


def main():
    ap = argparse.ArgumentParser(description="由 Mesh3d.dat + flow3d.dat 估算 CFL 数")
    ap.add_argument("--dir", default=".", help="算例目录 (含 control.ec, Mesh3d.dat)")
    ap.add_argument("--mesh", default="Mesh3d.dat")
    ap.add_argument("--flow", default="flow3d.dat")
    ap.add_argument("--dt", type=float, default=None, help="时间步长 (默认取 control.ec 中的 dt_global)")
    ap.add_argument("--freestream", action="store_true",
                    help="不读 flow3d.dat, 用均匀来流 (rho=1, |V|=1, T=1) 估算")
    ap.add_argument("--top", type=int, default=10, help="列出 CFL 最大的块数")
    args = ap.parse_args()

    os.chdir(args.dir)
    nml = read_namelist("control.ec")
    # 默认值与 sub_read_parameter.f90 一致
    para = dict(ma=1.0, re=1000.0, gamma=1.4, prl=0.7, t_inf=288.15, if_viscous=1,
                iflag_local_dt=0, cfl=1.0, dt_global=1.0e-4, dtmax=10.0, dtmin=1.0e-9,
                time_method=0)
    for key in para:
        if key in nml:
            para[key] = type(para[key])(float(nml[key]))
    local_dt = para["iflag_local_dt"] == 1
    dt = args.dt if args.dt is not None else para["dt_global"]

    print(f"Ma={para['ma']}  Re={para['re']}  If_viscous={para['if_viscous']}  "
          f"Time_Method={para['time_method']}  Iflag_local_dt={para['iflag_local_dt']}")
    if local_dt and args.dt is None:
        print(f"当地时间步长: dt = CFL({para['cfl']})*Vol/L, 并受 dtmax={para['dtmax']}, dtmin={para['dtmin']} 截断;"
              " 下面给出截断后的实际 CFL")
    else:
        print(f"全局时间步长 dt = {dt}")

    fm = FortranFile(args.mesh)
    ff = None if args.freestream else FortranFile(args.flow)
    nb = int(fm.read("<i4")[0])
    dims = fm.read("<i4").reshape(nb, 3)

    all_cfl, all_cflc, all_vol = [], [], []
    blk_info = []
    for m in range(nb):
        nx, ny, nz = (int(a) for a in dims[m])
        xyz = fm.read("<f8").reshape(nx, ny, nz, 3, order="F")
        if ff is None:
            one = np.ones((nx - 1, ny - 1, nz - 1))
            rho, u, v, w, T = one, one, 0 * one, 0 * one, one  # |V|=1, 方向不影响量级
        else:
            U = ff.read("<f8").reshape(nx + 1, ny + 1, nz + 1, 5, order="F")
            U = U[1:nx, 1:ny, 1:nz, :]  # 去掉 Ghost Cell
            rho, u, v, w, T = (U[..., i] for i in range(5))

        geo = block_geometry(xyz[..., 0], xyz[..., 1], xyz[..., 2])
        vol = geo[3]
        Lc, Lv = spectral_radius(geo, vol, rho, u, v, w, T, para)

        if local_dt and args.dt is None:
            dtl = np.clip(para["cfl"] * vol / (Lc + Lv), para["dtmin"], para["dtmax"])
        else:
            dtl = dt
        cfl = dtl * (Lc + Lv) / vol
        cflc = dtl * Lc / vol

        idx = np.unravel_index(np.argmax(cfl), cfl.shape)
        P = xyz[idx[0]:idx[0] + 2, idx[1]:idx[1] + 2, idx[2]:idx[2] + 2].reshape(-1, 3).mean(axis=0)
        blk_info.append((cfl[idx], cflc[idx], m + 1, tuple(int(i) + 1 for i in idx), P, np.min(vol)))
        all_cfl.append(cfl.ravel())
        all_cflc.append(cflc.ravel())
        all_vol.append(vol.ravel())
        print(f"\r  block {m + 1}/{nb}", end="", file=sys.stderr, flush=True)
    print(file=sys.stderr)
    fm.close()
    if ff is not None:
        ff.close()

    cfl = np.concatenate(all_cfl)
    cflc = np.concatenate(all_cflc)
    vol = np.concatenate(all_vol)
    if np.any(vol <= 0):
        print(f"警告: {np.sum(vol <= 0)} 个单元体积 <= 0")

    pct = [50, 90, 99, 99.9]
    print(f"\n单元总数: {cfl.size}")
    print(f"{'':10s}{'max':>12s}{'min':>12s}{'mean':>12s}" + "".join(f"{'P' + str(p):>12s}" for p in pct))
    for name, a in (("CFL", cfl), ("CFL_conv", cflc)):
        print(f"{name:10s}{a.max():12.4g}{a.min():12.4g}{a.mean():12.4g}"
              + "".join(f"{np.percentile(a, p):12.4g}" for p in pct))
    for lim in (1.0, 10.0, 100.0):
        print(f"  CFL > {lim:<6g}: {np.sum(cfl > lim) / cfl.size * 100:8.4f} %")

    blk_info.sort(key=lambda t: -t[0])
    print(f"\nCFL 最大的 {min(args.top, nb)} 个块 (i,j,k 为单元编号, 从 1 开始):")
    print(f"{'block':>6s}{'CFL_max':>12s}{'CFL_conv':>12s}{'(i,j,k)':>18s}{'x,y,z':>40s}{'min Vol':>12s}")
    for c, cc_, b, ijk, P, vmin in blk_info[:args.top]:
        print(f"{b:6d}{c:12.4g}{cc_:12.4g}{str(ijk):>18s}"
              f"{'(%.4g, %.4g, %.4g)' % tuple(P):>40s}{vmin:12.4g}")


if __name__ == "__main__":
    main()
