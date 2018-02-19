#!/usr/bin/env python

import os
import numpy as np
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import subprocess

NMODEL = 228
DT = 1.0  # write interval in [timeunit]
NT = 100

ONEDAY = 8.64  # [timeunit/day] a46p51

def main():
    np.set_printoptions(formatter={'float': '{: 0.5f}'.format})
    mkdirs()
    fname = "evol_field_tlm.dat"
    all_blv = np.empty((NT, NMODEL, NMODEL))
    blv = np.random.normal(0.0, 1.0, (NMODEL, NMODEL))
    blv, ble = orth_norm_vectors(blv)
    all_ble = np.zeros((NT, NMODEL))
    for i in range(0, NT):
        traj, tlm = read_file_part(fname, i)
        blv = np.dot(tlm[:, :], blv)
        blv, ble = orth_norm_vectors(blv)
        all_ble[i, :] = ble[:] / (DT / ONEDAY)
        all_blv[i, :, :] = blv[:, :]
    plt.plot(np.mean(all_ble, axis=0))
    print(np.mean(all_ble, axis=0))
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

def read_file_part(fname, it):
    dbyte = it * (NMODEL ** 2 + NMODEL) * 8
    x = np.memmap(fname, offset=dbyte, dtype=np.float64, mode='r', shape=(NMODEL ** 2 + NMODEL))
    traj = x[:NMODEL]
    tlm = x[NMODEL:].reshape((NMODEL, NMODEL))
    return traj, tlm

if __name__ == "__main__":
    main()
