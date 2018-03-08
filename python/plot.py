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

def get_grid_val(waves, x, y, elem):
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
            gridval *= L ** 2 * f0
        return gridval

    assert waves.__class__ == np.ndarray
    assert waves.shape == (36,)
    assert y.__class__ in [float, np.float32, np.float64]
    assert x.__class__ in [float, np.float32, np.float64]

    n = 1.5
    na = 10
    no = 8
    R = 287.0
    L = 5000000.0 / np.pi

    if elem == "a_psi":
        return get_atm(False)
    elif elem == "a_tmp":
        return get_atm(True)
    elif elem == "o_psi":
        return get_ocn(False)
    elif elem == "o_tmp":
        return get_ocn(True)
    else:
        raise Exception("get_grid_val overflow. Element %s not found." % elem)

def test_get_grid_val():
    n = 1.5
    state = model_state_exsample()
    x = 1.2 * np.pi / n  # 0.0 <= x <= 2.0 * pi / n
    y = 0.8 * np.pi      # 0.0 <= y <= pi

    # get_grid_val() returns one of four variables
    # {atmosphere|ocean} x {streamfunction|temperature} at point (x, y)
    # unit: [m^2/s] for streamfunction and [K] for temperature
    a_psi = get_grid_val(state, x, y, "a_psi")
    a_tmp = get_grid_val(state, x, y, "a_tmp")
    o_psi = get_grid_val(state, x, y, "o_psi")
    o_tmp = get_grid_val(state, x, y, "o_tmp")
    print(a_psi, a_tmp, o_psi, o_tmp)

def all_reconstruct_grid(waves, nx, ny):
    n = 1.5
    x_grid = np.linspace(0, 2.0 * np.pi / n, nx)
    y_grid = np.linspace(np.pi, 0, ny)
    a_psi = np.empty((ny, nx))
    a_tmp = np.empty((ny, nx))
    o_psi = np.empty((ny, nx))
    o_tmp = np.empty((ny, nx))
    for ix in range(nx):
        for iy in range(ny):
            a_psi[iy, ix] = get_grid_val(waves, x_grid[ix], y_grid[iy], "a_psi")
            a_tmp[iy, ix] = get_grid_val(waves, x_grid[ix], y_grid[iy], "a_tmp")
            o_psi[iy, ix] = get_grid_val(waves, x_grid[ix], y_grid[iy], "o_psi")
            o_tmp[iy, ix] = get_grid_val(waves, x_grid[ix], y_grid[iy], "o_tmp")
    o_psi -= np.mean(o_psi)
    return a_psi, a_tmp, o_psi, o_tmp

def model_state_exsample():
    xini = np.array([
        4.695340259215241E-002,
        2.795833230987369E-002,
        -2.471191763590483E-002,
        -7.877635082773315E-003,
        -4.448292568544942E-003,
        -2.756238610924190E-002,
        -4.224400051368891E-003,
        5.914241112882518E-003,
        -1.779437742222920E-004,
        5.224450720394076E-003,
        4.697982667229096E-002,
        5.149282577209392E-003,
        -1.949084549066326E-002,
        4.224006062949761E-004,
        -1.247786759371923E-002,
        -9.825952138046594E-003,
        -2.610941795170075E-005,
        2.239286581216401E-003,
        -7.891896725509534E-004,
        7.470171905055880E-004,
        -9.315932162526787E-007,
        3.650179005106874E-005,
        1.064122403269511E-006,
        3.937836448211443E-008,
        -2.208288760403859E-007,
        -3.753762121228048E-006,
        -7.105126469908465E-006,
        1.518110190916469E-008,
        -5.773178576933025E-004,
        0.187369278208256,
        1.369868543156558E-003,
        7.023608700166264E-002,
        -4.539810680860224E-004,
        -1.882650440363933E-003,
        -3.900412687995408E-005,
        -1.753655087903711E-007])
    return xini

if __name__ == "__main__":
    test_get_grid_val()
