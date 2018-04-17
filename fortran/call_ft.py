#!/usr/bin/env python
#vim:fileencoding=utf8
from ctypes import *
import numpy as np

n = 36

add_np = np.ctypeslib.load_library("step_maooam.so", ".")
add_np.step_maooam_.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64)]
add_np.step_maooam_.restype = c_void_p

x0 = np.random.randn(n) * 0.01
print(x0)
dt = np.array([0.01])
add_np.step_maooam_(x0, dt)
print(x0)

