-- ao_mpi.lua
-- (C) 2014-2016 Lesley De Cruz & Jonathan Demaeyer
-- See LICENSE.txt for license information.

------------------------------------------------------------------------
-- This file implements ao_mpi, a function which returns a parallel version of
-- ao, the function which computes the model tendencies. 
--
-- The function returned by ao_mpi is a rank-specific function which
-- distributes the calculation of the model tendencies over multiple cores.
--
-- It uses the lua-mpi module by Peter Colberg: https://colberg.org/lua-mpi/
------------------------------------------------------------------------

local mpi = require("mpi")
local tensor = require("tensor")

local inprod = require("inprod_analytic")
local n = inprod.natm*2+inprod.noc*2 + 1

-- Tensor with the coefficients of the nonlinear (polynomial) system of
-- diffeqs. N.B.: the model tensor is fixed, so changes to model parameters
-- will only be effective up to this point!
local aotensor = tensor.simplify_coo(require("aotensor"))

--- Get a rank-specific function which calculates part of the tensor contraction
-- and reduces the result to all cores.
local function ao_mpi(rank,size,comm)
  -- partition the tensor into subtensors for each core.
  local nelem_per_core = math.ceil(#aotensor/(size))
  local rank_offset = nelem_per_core*rank 
  local rank_aotensor = {}
  for i=1,nelem_per_core do
    -- Last core will generally not be filled up
    rank_aotensor[i] = aotensor[rank_offset+i]
  end
  rank_aotensor = tensor.coo_to_fficoo(rank_aotensor)

  local sparse_mul3 = tensor.sparse_mul3
  local sum, allreduce, in_place, double = mpi.sum, mpi.allreduce, mpi.in_place, mpi.double

  --- Rank-specific function that calculates the time derivative of the n variables.
  -- The function reduces to a sparse tensor contraction due to the bilinear
  -- nature of the equations.
  -- @function ao
  -- @param t time
  -- @param y array with variables at time t
  -- @param buf n-array (buffer) to store derivatives.
  return function(t,y,buf)
    buf = sparse_mul3(rank_aotensor,y,y,buf)
    allreduce(in_place, buf, n, double, sum, comm)
    return buf
  end
end

return ao_mpi
