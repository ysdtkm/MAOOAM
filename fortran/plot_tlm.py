#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import mydebug

N = 36
# N = 228   # atm 6x6 ocn 6x6
# N = 318   # atm 6x6 ocn 9x9
# N = 414   # atm 9x9 ocn 6x6
DT = 1.0  # write interval in [timeunit]
NT = 1000  # number of write. filesize = 8 * (N ** 2 + N) * NT [bytes]
ONEDAY = 8.64  # [timeunit/day] a46p51

def main():
    np.set_printoptions(formatter={'float': '{: 0.4f}'.format})
    mkdirs()
    fname = "/lustre/tyoshida/shrt/exec/m200/evol_field_tlm.dat"
    # all_blv = np.empty((NT, N, N))
    all_ble = np.zeros((NT, N))
    blv = np.random.normal(0.0, 1.0, (N, N))
    blv, ble = orth_norm_vectors(blv)
    tlms = np.empty((NT, N, N))
    for i in range(NT):
        traj, tlm = read_file_part(fname, i)
        tlms[i, :, :] = tlm
        blv = np.dot(tlm[:, :], blv)
        blv, ble = orth_norm_vectors(blv)
        all_ble[i, :] = ble[:] / (DT / ONEDAY)
        # all_blv[i, :, :] = blv[:, :]
    mydebug.dump_array(tlms)
    mean_ble = np.mean(all_ble[NT//2:, :], axis=0)
    plt.plot(mean_ble)
    print(mean_ble)
    plt.savefig("img/tmp.pdf")

def orth_norm_vectors(lv):
    # lv     <- np.array[N,N] : Lyapunov vectors (column)
    # return -> np.array[N,N] : orthonormalized LVs in descending order
    # return -> np.array[N]   : ordered Lyapunov Exponents
    q, r = np.linalg.qr(lv)
    le = np.zeros(N)
    eigvals = np.abs(np.diag(r))
    lv_orth = q
    le = np.log(eigvals)
    return lv_orth, le

def mkdirs():
    subprocess.run("rm -rf img", check=True, shell=True)
    subprocess.run("mkdir -p img", check=True, shell=True)

def read_file_part(fname, it):
    dbyte = it * (N ** 2 + N) * 8
    x = np.array(np.memmap(fname, offset=dbyte, dtype=np.float64, mode='r', shape=(N ** 2 + N)))
    assert x.shape == (N ** 2 + N, )
    traj = x[:N]
    tlm = x[N:].reshape((N, N))
    return traj, tlm

if __name__ == "__main__":
    main()
