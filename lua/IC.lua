-- IC.lua

-- Initial conditions for the extended model.

return {
-- psi
 0, -- psi(Nx=0, Ny=1, P=1, typ=A)
 0, -- psi(M=1, Nx=1, Ny=1, P=1, typ=K)
 0, -- psi(H=1, Nx=1, Ny=1, P=1, typ=L)
 0, -- psi(Nx=0, Ny=2, P=2, typ=A)
 0, -- psi(M=1, Nx=1, Ny=2, P=2, typ=K)
 0, -- psi(H=1, Nx=1, Ny=2, P=2, typ=L)
 0, -- psi(M=2, Nx=2, Ny=1, P=1, typ=K)
 0, -- psi(H=2, Nx=2, Ny=1, P=1, typ=L)
 0, -- psi(M=2, Nx=2, Ny=2, P=2, typ=K)
 0, -- psi(H=2, Nx=2, Ny=2, P=2, typ=L)
-- theta
 0, -- theta(Nx=0, Ny=1, P=1, typ=A)
 0, -- theta(M=1, Nx=1, Ny=1, P=1, typ=K)
 0, -- theta(H=1, Nx=1, Ny=1, P=1, typ=L)
 0, -- theta(Nx=0, Ny=2, P=2, typ=A)
 0, -- theta(M=1, Nx=1, Ny=2, P=2, typ=K)
 0, -- theta(H=1, Nx=1, Ny=2, P=2, typ=L)
 0, -- theta(M=2, Nx=2, Ny=1, P=1, typ=K)
 0, -- theta(H=2, Nx=2, Ny=1, P=1, typ=L)
 0, -- theta(M=2, Nx=2, Ny=2, P=2, typ=K)
 0, -- theta(H=2, Nx=2, Ny=2, P=2, typ=L)
-- A
 0, -- A(H=1, Nx=0.5, Ny=1, P=1)
 0, -- A(H=1, Nx=0.5, Ny=2, P=2)
 0, -- A(H=1, Nx=0.5, Ny=3, P=3)
 0, -- A(H=1, Nx=0.5, Ny=4, P=4)
 0, -- A(H=2, Nx=1, Ny=1, P=1)
 0, -- A(H=2, Nx=1, Ny=2, P=2)
 0, -- A(H=2, Nx=1, Ny=3, P=3)
 0, -- A(H=2, Nx=1, Ny=4, P=4)
-- T
 0, -- T(H=1, Nx=0.5, Ny=1, P=1)
 0, -- T(H=1, Nx=0.5, Ny=2, P=2)
 0, -- T(H=1, Nx=0.5, Ny=3, P=3)
 0, -- T(H=1, Nx=0.5, Ny=4, P=4)
 0, -- T(H=2, Nx=1, Ny=1, P=1)
 0, -- T(H=2, Nx=1, Ny=2, P=2)
 0, -- T(H=2, Nx=1, Ny=3, P=3)
 0, -- T(H=2, Nx=1, Ny=4, P=4)
}
