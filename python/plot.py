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
    nad, timd = read_file("evol_field.dat", 0.0)
    plot_time(nad, timd)
    plot_anime(nad)
    # plot_3d_trajectory(nad[:, 21], nad[:, 29], nad[:, 0])

def mkdirs():
    subprocess.run("rm -rf img", check=True, shell=True)
    subprocess.run("mkdir -p img", check=True, shell=True)
    subprocess.run("mkdir -p img/a_gph", check=True, shell=True)
    subprocess.run("mkdir -p img/a_t", check=True, shell=True)
    subprocess.run("mkdir -p img/o_psi", check=True, shell=True)
    subprocess.run("mkdir -p img/o_t", check=True, shell=True)

def read_file(file, discard):
    # return np.ndarray[time, NMODEL]
    with open(file, "r") as f:
        ar = f.read().split()
    ar2 = []
    n = len(ar)
    na = np.empty((n // (NMODEL + 1), NMODEL))
    tim = np.empty((n // (NMODEL + 1)))
    for i in range(n // (NMODEL + 1)):
        tim[i]   = ar[i * (NMODEL + 1)]
        na[i, :] = ar[i * (NMODEL + 1) + 1:(i + 1) * (NMODEL + 1)]
    nt = na.shape[0]
    nad = na[int(nt * discard):nt, :]
    timd = tim[int(nt * discard):nt]
    return nad, timd

def plot_time(nad, timd):
    for i in range(NMODEL):
        plt.plot(timd * 3.07e-4, nad[:, i])
        plt.xlabel("model year")
        plt.savefig("img/x_%02d.png" % i)
        plt.close()

def plot_snap(cmaxs, nad, i):
    it = i // PINTVL
    psia, ta, psio, to = all_reconstruct_grid(nad[i], 20, 20)
    datas = {"a_gph": psia * f0 / g, "a_t": ta, "o_psi": psio, "o_t": to}
    for cmp in cmaxs:
        title = "%s %04d" % (cmp, it)
        plot_matrix(datas[cmp], "img/%s/%s_%04d.png" % (cmp, cmp, it), title, cmaxs[cmp], ipol="none")

def plot_anime(nad):
    # cmaxs = {"a_gph": 500, "a_t": 20, "o_psi": 5e+5, "o_t": 40}  # DDV2016
    cmaxs = {"a_gph": 500, "a_t": 20, "o_psi": 3e+4, "o_t": 40}  # VL2016
    nt = nad.shape[0]
    for i in range(nt - PINTVL * ANIMAX, nt, PINTVL):
        plot_snap(cmaxs, nad, i)
    for dir in cmaxs:
        subprocess.run("convert -delay 8 -loop 0 ./img/%s/*.png ./img/%s/anime.gif" % (dir, dir), check=True, shell=True)

def plot_3d_trajectory(x, y, z):
    # 3D trajectory
    # plt.rcParams["font.size"] = 16
    fig = plt.figure()
    fig.subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.98, wspace=0.04, hspace=0.04)
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, label="traj", marker=".")
    # ax.set_xlim([-8e-4, 8e-4])
    # ax.set_ylim([-0.05, 0.2])
    # ax.set_zlim([0.02, 0.06])
    ax.set_xlabel("Psi o 2")
    ax.set_ylabel("Theta o 2")
    ax.set_zlabel("Psi a 1")
    plt.savefig("./img/traj.png")
    plt.close()

def plot_matrix(mat, out, title="", cmax=None, ipol="none"):
    # plt.rcParams["font.size"] = 14
    fig, ax = plt.subplots(1)
    if cmax is None:
        cmax = np.max(np.abs(mat))
    cm = ax.imshow(mat, cmap=plt.cm.RdBu_r, aspect=0.7, interpolation=ipol)
    cm.set_clim(-1.0 * cmax, cmax)
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    plt.colorbar(cm)
    plt.title(title)
    plt.savefig(out)
    plt.close()

def reconstruct_grid(waves, x, y, elem):
    # return np.arrays of atm_psi[K], atm_gph [m], ocn_psi[m2s-1], ocn_tmp [K]
    def get_atm(is_temp):
        types = ["A", "K", "L", "A", "K", "L", "K", "L", "K", "L"]
        hs = [0, 1, 1, 0, 1, 1, 2, 2, 2, 2]
        ps = [1, 1, 1, 2, 2, 2, 1, 1, 2, 2]
        gridval = 0.0
        for j in range(na):
            j_all = j + na if is_temp else j
            if types[j] == "A":
                gridval = gridval + waves[j_all] * np.sqrt(2.0) * np.cos(ps[j] * y)
            elif types[j] == "K":
                gridval = gridval + waves[j_all] * 2.0 * np.cos(hs[j] * n * x) * np.sin(ps[j] * y)
            else:
                gridval = gridval + waves[j_all] * 2.0 * np.sin(hs[j] * n * x) * np.sin(ps[j] * y)
        if is_temp:
            gridval *= (f0 ** 2 * L ** 2) / R
        else:
            gridval *= L ** 2 * f0
        return gridval

    def get_ocn(is_temp):
        hos = [1, 1, 1, 1, 2, 2, 2, 2]
        pos = [1, 2, 3, 4, 1, 2, 3, 4]
        gridval = 0.0
        for j in range(no):
            j_all = j + (na * 2 + no) if is_temp else j + na * 2
            gridval = gridval + waves[j_all] * 2.0 * np.sin(0.5 * hos[j] * n * x) * np.sin(pos[j] * y)
        if is_temp:
            gridval *= (f0 ** 2 * L ** 2) / R
        else:
            gridval -= np.mean(gridval)
            gridval *= L ** 2 * f0
        return gridval

    assert waves.__class__ == np.ndarray
    assert y.__class__ in [float, np.float32, np.float64]
    assert x.__class__ in [float, np.float32, np.float64]
    assert waves.shape == (36,)

    n = 1.5
    na = 10
    no = 8
    R = 287.0
    L = 5000000.0 / np.pi

    if elem == "a_gph":
        return get_atm(False)
    elif elem == "a_tmp":
        return get_atm(True)
    elif elem == "o_psi":
        return get_ocn(False)
    elif elem == "o_tmp":
        return get_ocn(True)
    else:
        raise Exception("reconstruct_grid overflow. Element %s not found." % elem)

def all_reconstruct_grid(waves, nx, ny):
    n = 1.5
    x_grid = np.linspace(0, 2.0 * np.pi / n, nx)
    y_grid = np.linspace(np.pi, 0, ny)
    a_gph = np.empty((ny, nx))
    a_tmp = np.empty((ny, nx))
    o_psi = np.empty((ny, nx))
    o_tmp = np.empty((ny, nx))
    for ix in range(nx):
        for iy in range(ny):
            a_gph[iy, ix] = reconstruct_grid(waves, x_grid[ix], y_grid[iy], "a_gph")
            a_tmp[iy, ix] = reconstruct_grid(waves, x_grid[ix], y_grid[iy], "a_tmp")
            o_psi[iy, ix] = reconstruct_grid(waves, x_grid[ix], y_grid[iy], "o_psi")
            o_tmp[iy, ix] = reconstruct_grid(waves, x_grid[ix], y_grid[iy], "o_tmp")
    return a_gph, a_tmp, o_psi, o_tmp

if __name__ == "__main__":
    main()
