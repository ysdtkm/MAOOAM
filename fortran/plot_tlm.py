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
oneday = 8.64  # [timeunit/day] a46p51

def main():
    mkdirs()
    traj, tlm, tim = read_file("evol_field_tlm.dat")

    nt = len(tim)
    dt = tim[1] - tim[0]

    all_blv = np.empty((nt, NMODEL, NMODEL))
    blv = np.random.normal(0.0, 1.0, (NMODEL, NMODEL))
    blv, ble = orth_norm_vectors(blv)
    all_ble = np.zeros((nt, NMODEL))
    for i in range(0, nt):
        blv = np.dot(tlm[i, :, :], blv)
        blv, ble = orth_norm_vectors(blv)
        all_ble[i, :] = ble[:] / (dt / oneday)
        all_blv[i, :, :] = blv[:, :]
    plt.plot(np.mean(all_ble, axis=0))
    plt.savefig("img/tmp.pdf")

def orth_norm_vectors(lv):
  # lv     <- np.array[NMODEL,NMODEL] : Lyapunov vectors (column)
  # return -> np.array[NMODEL,NMODEL] : orthonormalized LVs in descending order
  # return -> np.array[NMODEL]      : ordered Lyapunov Exponents

  q, r = np.linalg.qr(lv)

  le = np.zeros(NMODEL)
  eigvals = np.abs(np.diag(r))

  # for t-continuity, align
  for i in range(NMODEL):
    inner_prod = np.dot(q[:,i].T, lv[:,i])
    if (inner_prod < 0):
      q[:,i] *= -1.0

  lv_orth = q
  le = np.log(eigvals)
  return lv_orth, le

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
    nrec = NMODEL ** 2 + NMODEL + 1

    na = np.empty((n // nrec, NMODEL))
    ntl = np.empty((n // nrec, NMODEL ** 2))
    tim = np.empty((n // nrec))

    for i in range(n // nrec):
        tim[i]   = ar[i * nrec]
        na[i, :] = ar[i * nrec + 1:i * nrec + NMODEL + 1]
        ntl[i, :] = ar[i * nrec + NMODEL + 1:i * nrec + NMODEL ** 2 + NMODEL + 1]
    ntl = ntl.reshape((n // nrec, NMODEL, NMODEL))
    return na, ntl, tim

if __name__ == "__main__":
    main()
