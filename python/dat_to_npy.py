#!/usr/bin/env python

import numpy as np

N = 37

def read_file(fname_from, fname_to):
    # return np.ndarray[time, NMODEL]

    with open(fname_from, "r") as f:
        ar = f.read().split()

    ar2 = []
    n = len(ar)
    na = np.empty((n // N, N))
    for i in range(n // N):
        na[i, :] = ar[i * N:(i + 1) * N]

    np.save(fname_to, na)

read_file("evol_field.dat", "evol_field.npy")
