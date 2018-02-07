#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import subprocess
from mpl_toolkits.mplot3d import Axes3D

NMODEL = 36
PINTVL = 10
ANIMAX = 3

f0 = 1.032e-4
g = 9.81

def main():
    mkdirs()
    nad, timd = read_file("evol_field_tlm.dat")
    print(nad.shape)
    # plot_time(nad, timd)
    # plot_anime(nad)
    # plot_3d_trajectory(nad[:, 21], nad[:, 29], nad[:, 0])

def mkdirs():
    subprocess.run("rm -rf img", check=True, shell=True)
    subprocess.run("mkdir -p img", check=True, shell=True)
    subprocess.run("mkdir -p img/a_gph", check=True, shell=True)
    subprocess.run("mkdir -p img/a_t", check=True, shell=True)
    subprocess.run("mkdir -p img/o_psi", check=True, shell=True)
    subprocess.run("mkdir -p img/o_t", check=True, shell=True)

def read_file(file):
    # return np.ndarray[time, NMODEL]
    with open(file, "r") as f:
        ar = f.read().split()
    ar2 = []
    n = len(ar)
    na = np.empty((n // (NMODEL ** 2 + NMODEL + 1), NMODEL ** 2 + NMODEL))
    tim = np.empty((n // (NMODEL ** 2 + NMODEL + 1)))
    for i in range(n // (NMODEL ** 2 + NMODEL + 1)):
        tim[i]   = ar[i * (NMODEL ** 2 + NMODEL + 1)]
        na[i, :] = ar[i * (NMODEL ** 2 + NMODEL + 1) + 1:(i + 1) * (NMODEL ** 2 + NMODEL + 1)]
    nt = na.shape[0]
    return na, tim

# def plot_time(nad, timd):
#     for i in range(NMODEL):
#         plt.plot(timd * 3.07e-4, nad[:, i])
#         plt.xlabel("model year")
#         plt.savefig("img/x_%02d.png" % i)
#         plt.close()
# 
# def plot_snap(cmaxs, nad, i):
#     it = i // PINTVL
#     psia, ta, psio, to = reconstruct_grid(nad[i], 20, 20)
#     datas = {"a_gph": psia * f0 / g, "a_t": ta, "o_psi": psio, "o_t": to}
#     for cmp in cmaxs:
#         title = "%s %04d" % (cmp, it)
#         plot_matrix(datas[cmp], "img/%s/%s_%04d.png" % (cmp, cmp, it), title, cmaxs[cmp], ipol="none")
# 
# def plot_anime(nad):
#     cmaxs = {"a_gph": 500, "a_t": 20, "o_psi": 5e+5, "o_t": 40}
#     nt = nad.shape[0]
#     for i in range(nt - PINTVL * ANIMAX, nt, PINTVL):
#         plot_snap(cmaxs, nad, i)
#     for dir in cmaxs:
#         subprocess.run("convert -delay 8 -loop 0 ./img/%s/*.png ./img/%s/anime.gif" % (dir, dir), check=True, shell=True)
# 
# def plot_3d_trajectory(x, y, z):
#     # 3D trajectory
#     # plt.rcParams["font.size"] = 16
#     fig = plt.figure()
#     fig.subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.98, wspace=0.04, hspace=0.04)
#     ax = fig.add_subplot(111, projection='3d')
#     ax.scatter(x, y, z, label="traj", marker=".")
#     ax.set_xlim([-8e-4, 8e-4])
#     ax.set_ylim([-0.05, 0.2])
#     ax.set_zlim([0.02, 0.06])
#     ax.set_xlabel("Psi o 2")
#     ax.set_ylabel("Theta o 2")
#     ax.set_zlabel("Psi a 1")
#     plt.savefig("./img/traj.png")
#     plt.close()
# 
# def plot_matrix(mat, out, title="", cmax=None, ipol="none"):
#     # plt.rcParams["font.size"] = 14
#     fig, ax = plt.subplots(1)
#     if cmax is None:
#         cmax = np.max(np.abs(mat))
#     cm = ax.imshow(mat, cmap=plt.cm.RdBu_r, aspect=0.7, interpolation=ipol)
#     cm.set_clim(-1.0 * cmax, cmax)
#     x0, x1 = ax.get_xlim()
#     y0, y1 = ax.get_ylim()
#     plt.colorbar(cm)
#     plt.title(title)
#     plt.savefig(out)
#     plt.close()
# 
# def reconstruct_grid(waves, nx, ny):
#     def fa(p):
#         # return np.ndarray[ny, nx]
#         return np.sqrt(2) * np.cos(p * y_grid)
# 
#     def fk(m, p):
#         # return np.ndarray[ny, nx]
#         return 2.0 * np.cos(m * n * x_grid) * np.sin(p * y_grid)
# 
#     def fl(h, p):
#         # return np.ndarray[ny, nx]
#         return 2.0 * np.sin(h * n * x_grid) * np.sin(p * y_grid)
# 
#     def phio(ho, po):
#         # return np.ndarray[ny, nx]
#         return 2.0 * np.sin(0.5 * ho * n * x_grid) * np.sin(po * y_grid)
# 
#     def atm(is_t):
#         # return np.ndarray[ny, nx]
#         gridval = 0.0
#         for j in range(na):
#             j_all = j + na if is_t else j
#             if types[j] == "A":
#                 gridval = gridval + waves[j_all] * fa(ps[j])
#             elif types[j] == "K":
#                 gridval = gridval + waves[j_all] * fk(hs[j], ps[j])
#             else:
#                 gridval = gridval + waves[j_all] * fl(hs[j], ps[j])
#         if is_t:
#             # gridval *= (2.0 * f0 / R)
#             gridval *= (f0 ** 2 * L ** 2) / R
#         else:
#             gridval *= L ** 2 * f0
#         return gridval
# 
#     def ocn(is_t):
#         # return np.ndarray[ny, nx]
#         gridval = 0.0
#         for j in range(no):
#             j_all = j + (na * 2 + no) if is_t else j + na * 2
#             gridval = gridval + waves[j_all] * phio(hos[j], pos[j])
#         if is_t:
#             gridval *= (f0 ** 2 * L ** 2) / R
#         else:
#             gridval -= np.mean(gridval)
#             gridval *= L ** 2 * f0
#         return gridval
# 
#     n = 1.5
#     x_grid = np.empty((ny, nx))
#     x_grid[:, :] = np.linspace(0, 2.0 * np.pi / n, nx)[np.newaxis, :]
#     y_grid = np.empty((ny, nx))
#     y_grid[:, :] = np.linspace(np.pi, 0, ny)[:, np.newaxis]
# 
#     na = 10
#     no = 8
#     R = 287
#     L = 5000000 / np.pi
#     types = ["A", "K", "L", "A", "K", "L", "K", "L", "K", "L"]
#     hs = [0, 1, 1, 0, 1, 1, 2, 2, 2, 2]
#     ps = [1, 1, 1, 2, 2, 2, 1, 1, 2, 2]
#     hos = [1, 1, 1, 1, 2, 2, 2, 2]
#     pos = [1, 2, 3, 4, 1, 2, 3, 4]
# 
#     return atm(False), atm(True), ocn(False), ocn(True)

if __name__ == "__main__":
    main()
