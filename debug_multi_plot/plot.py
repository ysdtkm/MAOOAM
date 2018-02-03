#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
mpl.use('PDF')
import matplotlib.pyplot as plt
import subprocess
from mpl_toolkits.mplot3d import Axes3D
import multiprocessing, functools

NMODEL = 36
PINTVL = 10
ANIMAX = 10

f0 = 1.032e-4
g = 9.81

def main():
    mkdirs()
    nad, timd = read_npy("evol_field.npy", 0.0, 0.1)
    # plot_time(nad, timd)
    plot_anime(nad)

def mkdirs():
    subprocess.run("rm -rf img", check=True, shell=True)
    subprocess.run("mkdir -p img", check=True, shell=True)
    subprocess.run("mkdir -p img/a_gph", check=True, shell=True)
    subprocess.run("mkdir -p img/a_t", check=True, shell=True)
    subprocess.run("mkdir -p img/o_psi", check=True, shell=True)
    subprocess.run("mkdir -p img/o_t", check=True, shell=True)

def read_npy(file, st, ed):
    na = np.load(file)
    nt = na.shape[0]
    nad = na[int(nt * st):int(nt * ed), 1:]
    timd = na[int(nt * st):int(nt * ed), 0]
    return nad, timd

def plot_time(nad, timd):
    for i in range(NMODEL):
        plt.plot(timd * 3.07e-4, nad[:, i])
        plt.xlabel("model year")
        plt.savefig("img/x_%02d.png" % i)
        plt.close()

def plot_snap(cmaxs, nad, i):
    it = i // PINTVL
    psia, ta, psio, to = reconstruct_grid(nad[i], 20, 20)
    datas = {"a_gph": psia * f0 / g, "a_t": ta, "o_psi": psio, "o_t": to}
    for cmp in cmaxs:
        title = "%s %04d" % (cmp, it)
        plot_matrix(datas[cmp], "img/%s/%s_%04d.png" % (cmp, cmp, it), title, cmaxs[cmp], ipol="none")

def plot_anime(nad):
    cmaxs = {"a_gph": 500, "a_t": 20, "o_psi": 5e+5, "o_t": 40}
    nt = nad.shape[0]
    p = multiprocessing.Pool(4)
    p.map(functools.partial(plot_snap, cmaxs, nad), range(nt - PINTVL * ANIMAX, nt, PINTVL))
    p.close()
    for dir in cmaxs:
        subprocess.run("convert -delay 15 -loop 0 ./img/*/*.png ./img/anime.gif", check=True, shell=True)

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

def reconstruct_grid(waves, nx, ny):
    def fa(p):
        # return np.ndarray[ny, nx]
        return np.sqrt(2) * np.cos(p * y_grid)

    def fk(m, p):
        # return np.ndarray[ny, nx]
        return 2.0 * np.cos(m * n * x_grid) * np.sin(p * y_grid)

    def fl(h, p):
        # return np.ndarray[ny, nx]
        return 2.0 * np.sin(h * n * x_grid) * np.sin(p * y_grid)

    def phio(ho, po):
        # return np.ndarray[ny, nx]
        return 2.0 * np.sin(0.5 * ho * n * x_grid) * np.sin(po * y_grid)

    def atm(is_t):
        # return np.ndarray[ny, nx]
        gridval = 0.0
        for j in range(na):
            j_all = j + na if is_t else j
            if types[j] == "A":
                gridval = gridval + waves[j_all] * fa(ps[j])
            elif types[j] == "K":
                gridval = gridval + waves[j_all] * fk(hs[j], ps[j])
            else:
                gridval = gridval + waves[j_all] * fl(hs[j], ps[j])
        if is_t:
            # gridval *= (2.0 * f0 / R)
            gridval *= (f0 ** 2 * L ** 2) / R
        else:
            gridval *= L ** 2 * f0
        return gridval

    def ocn(is_t):
        # return np.ndarray[ny, nx]
        gridval = 0.0
        for j in range(no):
            j_all = j + (na * 2 + no) if is_t else j + na * 2
            gridval = gridval + waves[j_all] * phio(hos[j], pos[j])
        if is_t:
            gridval *= (f0 ** 2 * L ** 2) / R
        else:
            gridval -= np.mean(gridval)
            gridval *= L ** 2 * f0
        return gridval

    n = 1.5
    x_grid = np.empty((ny, nx))
    x_grid[:, :] = np.linspace(0, 2.0 * np.pi / n, nx)[np.newaxis, :]
    y_grid = np.empty((ny, nx))
    y_grid[:, :] = np.linspace(np.pi, 0, ny)[:, np.newaxis]

    na = 10
    no = 8
    R = 287
    L = 5000000 / np.pi
    types = ["A", "K", "L", "A", "K", "L", "K", "L", "K", "L"]
    hs = [0, 1, 1, 0, 1, 1, 2, 2, 2, 2]
    ps = [1, 1, 1, 2, 2, 2, 1, 1, 2, 2]
    hos = [1, 1, 1, 1, 2, 2, 2, 2]
    pos = [1, 2, 3, 4, 1, 2, 3, 4]

    return atm(False), atm(True), ocn(False), ocn(True)

if __name__ == "__main__":
    main()
