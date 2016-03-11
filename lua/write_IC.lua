-- write_IC.lua
-- (C) 2014-2016 Lesley De Cruz & Jonathan Demaeyer
-- See LICENSE.txt for license information.

------------------------------------------------------------------------
-- The function write_IC writes an initial condition template in IC.lua if the
-- file does not yet exist. It then also updates the file called
-- translation.txt, which contains a translation of the linear array of
-- variables into the specific ocean and atmosphere variable names.
------------------------------------------------------------------------

local function tabstr(t)
  local str={}
  for k,v in pairs(t) do str[#str+1] = k.."="..v end
  table.sort(str)
  return str
end

local function write_IC(inprod, initvec)
  -- write IC template if necessary, return true if it is written.
  local IC = io.open("IC.lua")
  if IC then return false end
  IC=io.open("IC.lua","w")
  assert(IC,"IC.lua is not present and not writeable.")
  IC:write("-- IC.lua\n\n-- Initial conditions for the extended model.\n\nreturn {\n")
  initvec = initvec or setmetatable({},{__index=function() return 0 end})
  local i_vec = 1
  for _,var in ipairs{"psi","theta"} do
    IC:write("-- ",var,"\n")
    for _,wavenums in ipairs(inprod.get_awavenum) do
      IC:write(" "..initvec[i_vec].. ", -- ",var,"(",table.concat(tabstr(wavenums),", "),")\n")
      i_vec = i_vec + 1
    end
  end
  for _,var in ipairs{"A","T"} do
    IC:write("-- ",var,"\n")
    for _,wavenums in ipairs(inprod.get_owavenum) do
      IC:write(" "..initvec[i_vec]..", -- ",var,"(",table.concat(tabstr(wavenums),", "),")\n")
      i_vec = i_vec + 1
    end
  end
  IC:write("}\n")
  IC:close()

  local tab=io.open("translation.txt","w")  -- always update translation.txt
  assert(IC,"Translation file is not present and not writeable.")
  tab:write("-- Translation table for the extended model :\n\n")
  for _,var in ipairs{"psi","theta"} do
    tab:write("-- ",var,"\n")
    local i=0
    for _,wavenums in ipairs(inprod.get_awavenum) do
      i=i+1
      tab:write(var,"_",string.format("%d",i),"  <==>  ",table.concat(tabstr(wavenums),", "),"\n")
    end
    tab:write("\n")
  end
  for _,var in ipairs{"A","T"} do
    tab:write("-- ",var,"\n")
    local i=0
    for _,wavenums in ipairs(inprod.get_owavenum) do
      i=i+1
      tab:write(var,"_",string.format("%d",i),"  <==>  ",table.concat(tabstr(wavenums),", "),"\n")
    end
    tab:write("\n")
  end
  tab:write("\n")
  tab:close()
  return true
end

return write_IC
