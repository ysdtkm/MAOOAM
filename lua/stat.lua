-- stat.lua
-- (C) 2013 Lesley De Cruz
-- See LICENSE.txt for license information.

------------------------------------------------------------------------
-- Stable statistics accumulator for n-arrays.
------------------------------------------------------------------------

--- Create a statistics accumulator for n-arrays.
-- @param n length of the array
-- @return statistics object for n-arrays
return function(n)
  local array = require("array")(n)
  local m, mprev, v, i = array(), array(), array(), 0
  local delta, delta2 = array(), array()
  --- add the item to the incremental mean and variance and increase
  -- the iteration count
  -- @param item n-array
  local function acc(item)
    i = i + 1
    item:sub(m,delta)               -- delta = item - m
    m:add(delta:div(i,mprev),mprev) -- mprev  = m + delta/i
    mprev, m = m, mprev
    v:add(delta:mul(item:sub(m,delta2),delta),v) -- v = v + (item - mprev) * (item - m)
  end
  --- get the incremental mean value
  local function mean() return m end
  --- get the incremental variance (sigma^2)
  local function var()  return v/(i-1) end
  --- get the iteration count
  local function iter() return i end
  --- set the statistics and counter
  local function set(new_m, new_v, new_i)
    m, v, i = new_m, new_v*(new_i-1), new_i
  end
  --- reset counter and clear statistics
  local function reset()
    set(array(),array(),0)
  end
  --- @export
  return {acc=acc, mean=mean, var=var, iter=iter, reset=reset, set=set}
end
