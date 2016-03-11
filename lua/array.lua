-- array.lua
-- (C) 2013 Lesley De Cruz
-- See LICENSE.txt for license information.

------------------------------------------------------------------------
-- Arrays with operations.
--
-- This module provides a function which returns a specialized
-- constructor for ffi-based arrays of size n (zero-filled by default).
--
-- The arrays created with this constructor have overloaded (pointwise)
-- operators + (add), - (sub), * (mul), / (div) and unary - (unm), as
-- well as __tostring and #. The modulo operator % (mod) calculates the
-- scalar product.
--
-- However, as these mathematical operations are not allowed to modify
-- the original arrays, they create new arrays for each operation. If
-- you do not want this, use the functions add, sub, mul, smul, div, unm
-- provided as methods of an array object. These functions take an extra
-- argument, namely a result buffer, which in turn is returned by these
-- functions, so that the operations can be chained.
--
-- To clear (i.e. zero-fill) an array, the `clear()` method can be called.
--
-- e.g.
--     local narray = require("array")(10)
--     local x, y = narray(), narray()
--     x[2] = 29
--     y:add(x,y) -- add y to x
--     print(y)
--     y:clear()
------------------------------------------------------------------------

local ffi = require("ffi")
local dsz = ffi.sizeof("double")
local assert, format = assert, string.format
if not DEBUG then assert = function(v) return v end end
local cache = {}

--- Create a constructor for n-arrays.
-- @function array
-- @param n array length 
-- @return n-array constructor
-- @export
return function(n)
  if cache[n] then return cache[n] end
  --- Multiplication with a number. 
  -- @function nmul
  -- @param lhs n-array
  -- @param factor number
  -- @param res result n-array (buffer)
  -- @return result n-array
  local function nmul (lhs, factor, res)
    for i=0,n-1 do res.data[i] = lhs.data[i]*factor end
    return res
  end
  -- Array operators.
  -- Note that as these operators perform a single pass, it is safe to provide
  -- lhs or rhs as the result buffer.
  -- @export
  local array_ops = {
    --- Pointwise addition
    add = function(lhs, rhs, res)
      for i=0,n-1 do res.data[i] = lhs.data[i]+rhs.data[i] end
      return res
    end,
    --- Pointwise subtraction
    sub = function(lhs, rhs, res)
      for i=0,n-1 do res.data[i] = lhs.data[i]-rhs.data[i] end
      return res
    end,
    nmul = nmul,
    --- Pointwise product
    mul = function(lhs, rhs, res)
      if type(rhs)=="number" then
        return nmul(lhs, rhs, res)
      elseif type(lhs)=="number" then
        return nmul(rhs, lhs, res)
      else -- element-wise multiplication
        for i=0,n-1 do res.data[i] = lhs.data[i]*rhs.data[i] end
        return res
      end
    end,
    --- Scalar product
    smul = function(lhs, rhs)
      local sp = 0
      for i=0,n-1 do sp = sp + lhs.data[i]*rhs.data[i] end
      return sp
    end,
    --- Pointwise division
    div = function(lhs, rhs, res)
      if type(rhs)=="number" then
        for i=0,n-1 do res.data[i] = lhs.data[i]/rhs end
      elseif type(lhs)=="number" then
        for i=0,n-1 do res.data[i] = lhs/rhs.data[i] end
      else
        for i=0,n-1 do res.data[i] = lhs.data[i]/rhs.data[i] end
      end
      return res
    end,
    --- Unary minus
    unm = function(lhs, res)
      for i=0,n-1 do res.data[i] = -lhs.data[i] end
      return res
    end,
    --- Zero-fill the array
    clear = function(lhs)
      ffi.fill(lhs,n*dsz)
      return lhs
    end,
    -- Copy
    copy = function(target,source)
      ffi.copy(target.data,source.data,n*dsz)
      return target
    end,
  }
  local narray
  local vec_mt_n = {
    __len = function() return n end,
    __add = function(lhs, rhs)
      return array_ops.add(lhs, rhs, narray())
    end,
    __sub = function(lhs, rhs)
      return array_ops.sub(lhs, rhs, narray())
    end,
    __unm = function(lhs)
      return array_ops.unm(lhs,narray())
    end,
    -- mul and div are pointwise operators (% is scalar)
    __mul = function(lhs, rhs) -- pointwise product
      return array_ops.mul(lhs, rhs, narray())
    end,
    __div = function(lhs, rhs) -- pointwise division
      return array_ops.div(lhs, rhs, narray())
    end,
    __mod = function(lhs, rhs) -- scalar product
      return array_ops.smul(lhs, rhs)
    end,
    __index = function(s, key)
      return array_ops[key] or assert(key<n,key) and s.data[key]
    end,
    __newindex = function(s,key,val)
      assert(key<n,key)
      s.data[key]=val
    end,
    __tostring = function(y)
      local t = {}
      for i=0,n-1 do t[i+1] = format("%.12f",y[i]) end
      return table.concat(t,"\t")
    end
  }
  narray = ffi.typeof("struct {double data[$];}",n)
  ffi.metatype(narray,vec_mt_n)
  cache[n] = narray
  return narray
end
