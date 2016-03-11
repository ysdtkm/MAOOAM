-- Code to test the aotensor (against the fortran code)

local tt = require("tensor")
local out=io.open("test_aotensor.ref","w")
-- Simplify the coo (has been done in Fortran code as well)
local ao_coo = tt.simplify_coo(require("aotensor"))
-- sort it just like in the fortran code.
table.sort(ao_coo,function(a,b) return a[1] < b[1] end)

for _,coo in ipairs(ao_coo) do
  out:write("aotensor")
  for i=1,#coo-1 do
    out:write(string.format("[%d]",coo[i]))
  end
  out:write(string.format(" = % .5E\n",coo[#coo]))
end
out:close()

