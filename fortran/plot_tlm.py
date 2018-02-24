#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import subprocess

N = 36
# N = 228   # atm 6x6 ocn 6x6
# N = 318   # atm 6x6 ocn 9x9
# N = 414   # atm 9x9 ocn 6x6
DT = 1.0  # write interval in [timeunit]
NT = 1000  # number of write. filesize = 8 * (N ** 2 + N) * NT [bytes]
ONEDAY = 8.64  # [timeunit/day] a46p51
FNAME = "/lustre/tyoshida/shrt/exec/m200/evol_field_tlm.dat"
LOGNORM = True

def main():
    np.set_printoptions(formatter={'float': '{: 0.4f}'.format})
    mkdirs()
    trajs, gs, ms, rs = integ_forward_ginelli()
    cs = integ_backward_ginelli(rs)
    vs = obtain_clvs_ginelli(cs, gs)
    test_growth_rate(ms, vs)
    test_growth_rate_long(ms, vs, 100, [0, 1, 24, 29])

def integ_forward_ginelli():
    # trajs[i, :] and gs[i, :, :] are at time i
    # ms[i, :, :] and rs[i, :, :] are times i -> i + 1
    g = np.random.normal(0.0, 1.0, (N, N))
    g, r = np.linalg.qr(g)
    trajs = np.empty((NT, N))
    gs = np.empty((NT, N, N))
    ms = np.empty((NT, N, N))
    rs = np.empty((NT, N, N))
    for i in range(NT):
        traj, m = read_file_part(FNAME, i)
        trajs[i, :] = traj
        gs[i, :, :] = g
        g = np.dot(m, g)
        g, r = np.linalg.qr(g)
        ms[i, :, :] = m
        rs[i, :, :] = r
    return trajs, gs, ms, rs

def integ_backward_ginelli(rs):
    # rs[i, :, :] is times i -> i + 1
    # cs[i, :, :] is at time i
    cs = np.empty((NT, N, N))
    c = np.random.normal(0.0, 1.0, (N, N))
    dummy, c = np.linalg.qr(c)
    c = normalize_column(c)
    for i in reversed(range(NT)):
        r = rs[i, :, :]
        ri = np.linalg.inv(r)
        c = np.dot(ri, c)
        c = normalize_column(c)
        cs[i, :, :] = c
    return cs

def obtain_clvs_ginelli(cs, gs):
    # cs[i, :, :], gs[i, :, :], and vs[i, :, :] are at time i
    vs = np.empty((NT, N, N))
    for i in range(NT):
        c = cs[i, :, :]
        g = gs[i, :, :]
        v = np.dot(g, c)
        vs[i, :, :] = v
    return vs

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

def normalize_column(m):
    for i in range(N):
        c = m[:, i]
        lc = np.dot(c, c) ** 0.5
        m[:, i] /= lc
    return m

def test_growth_rate(ms, vs):
    rate = [0.0 for j in range(N)]
    for i in range(NT):
        v = vs[i, :, :]
        m = ms[i, :, :]
        vn = np.dot(m, v)
        for j in range(N):
            r = np.log(np.linalg.norm(vn[:, j]) / np.linalg.norm(v[:, j]))
            rate[j] += r / (DT / ONEDAY) / NT
    print("Growth rate:")
    print(rate)

def test_growth_rate_long(ms, vs, ntg, ilist):
    lengths = np.zeros((NT - ntg, ntg + 1, len(ilist)))
    for k, i in enumerate(ilist):
        for it in range(NT - ntg):
            vi = vs[it, :, i]
            lengths[it, 0, k] = np.linalg.norm(vi)
            for jt in range(ntg):
                vi = np.dot(ms[it + jt, :, :], vi[:])
                lengths[it, jt + 1, k] = np.linalg.norm(vi)
    growths = np.copy(lengths)
    for jt in range(ntg):
        growths[:, jt + 1, :] /= lengths[:, jt, :]
    growths[:, 0, :] = np.nan
    if LOGNORM:
        growth_mean = np.mean(np.log(growths), axis=0) / (DT / ONEDAY)
    else:
        growth_mean = np.log(np.mean(growths, axis=0)) / (DT / ONEDAY)
    plot_growth_rate(growth_mean, ilist)

def plot_growth_rate(growth_mean, ilist):
    # growth_mean.shape = (ntg + 1, len(ilist))
    for k, i in enumerate(ilist):
        plt.plot(growth_mean[:, k], label="CLV %02d" % (i + 1))
    plt.legend()
    plt.savefig("img/fig_8.png")
    plt.close()

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
