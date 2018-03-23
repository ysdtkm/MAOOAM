#!/usr/bin/env python3

import sys
import subprocess as sp
import numpy as np

def update_params_nml(coupl):
    assert coupl.__class__ in [int, float, np.float32, np.float64]
    fname = "/home/tak/prgm/MAOOAM/fortran/params.nml"
    line_num = 19
    with open(fname) as f:
        st = f.read().splitlines()
    st2 = []
    for i, l in enumerate(st):
        if i + 1 == line_num:
            assert "COUPLE = " in l
            l = "  COUPLE = %f" % coupl
        st2.append(l)
    st2 = "\n".join(st2)
    with open(fname, "w") as f:
        f.write(st2)

def loop():
    n = 11
    coupls = np.linspace(0.0, 1.0, n)
    sp.run("make")
    sp.run("mkdir -p dat", shell=True, check=True)
    for i in range(n):
        update_params_nml(coupls[i])
        sp.run("./run_true_and_tlm", shell=True, check=True)
        sp.run("mv evol_field_tlm.dat dat/evol_field_tlm_%.2f.dat" % coupls[i], shell=True, check=True)
    
if __name__ == "__main__":
    loop()
    
