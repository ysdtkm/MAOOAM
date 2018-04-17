#!/usr/bin/env python
from ctypes import *
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

n = 36

class Maooam_Fortran:
    dt = np.array([0.1])
    one = np.array([1.0])

    add_np = np.ctypeslib.load_library("step_maooam.so", ".")
    add_np.step_maooam_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64)]
    add_np.step_maooam_.restype = c_void_p

    def step(self, x0):
        y0 = np.concatenate((self.one, x0))
        self.add_np.step_maooam_(y0, self.dt)
        return y0[1:]

def main():
    nt = 10000
    mf = Maooam_Fortran()
    x0 = np.random.randn(n) * 0.01
    xhist = np.empty((nt, n))
    for i in range(10000):
        x0 = mf.step(x0)
        xhist[i, :] = x0
    plot_xhist(xhist)

def plot_xhist(data):
    cm = plt.imshow(data, aspect="auto")
    plt.colorbar(cm)
    plt.savefig("tmp.pdf")
    plt.close()

main()

