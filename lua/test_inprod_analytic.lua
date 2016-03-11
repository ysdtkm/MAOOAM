-- test_inprod_analytic.lua
-- (C) 2013 Lesley De Cruz
-- See LICENSE.txt for license information.

------------------------------------------------------------------------
-- This script writes out the inner products sparse tensor.
------------------------------------------------------------------------

local inprod = require("inprod_analytic")
local abs = math.abs
local real_eps = 2.2204460492503131e-16
local f = io.open("test_inprod_analytic.ref","w")
local write = function(...)
  f:write(...)
end

local function dump(vn,name)
  if tonumber(vn) and abs(vn)>real_eps then
    write(name," = ",string.format("% .5E\n",vn))
  end
  return
end

local natm, noc = inprod.natm, inprod.noc
for i=1,natm do
  for j=1,natm do
    for _,name in ipairs{"a","c"} do
      dump(inprod.atmos[name][i][j],
             name.."["..i.."]["..j.."]")
    end
    for k=1,natm do
      for _,name in ipairs{"b","g"} do
        dump(inprod.atmos[name][i][j][k],
               name.."["..i.."]["..j.."]["..k.."]")
      end
    end
  end
  for j=1,noc do
    for _,name in ipairs{"d","s"} do
      dump( inprod.atmos[name][i][j],
           name.."["..i.."]["..j.."]")
    end
  end
end

for i=1,noc do
  for j=1,noc do
    for _,name in ipairs{"M","N"} do
      dump(inprod.ocean[name][i][j],
      name.."["..i.."]["..j.."]")
    end
    for k=1,noc do
      for _,name in ipairs{"O","C"} do
        dump(inprod.ocean[name][i][j][k],
               name.."["..i.."]["..j.."]["..k.."]")
      end
    end
  end
  for j=1,natm do
    for _,name in ipairs{"K","W"} do
      dump(inprod.ocean[name][i][j],
      name.."["..i.."]["..j.."]")
    end
  end
end
