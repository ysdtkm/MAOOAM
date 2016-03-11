-- modeselection.lua
-- (C) 2014-2016 Lesley De Cruz & Jonathan Demaeyer
-- See LICENSE.txt for license information.

------------------------------------------------------------------------
-- This file must return a table with two fields: atmos and ocean, each of
-- which contain an array with the selected {Nx,Ny} wavenumbers of the modes
-- to be included in the respective components of MAOOAM.
-- e.g.
--    return {
--      atmos = {{1,1},{1,2},{2,1},{2,2}},
--      ocean = {{1,1},{1,2},{2,1},{2,2}},
--    }
-- 
-- It includes a utility function getmodes(Nx_max, Ny_max) which generates all
-- combinations of wavenumbers from {1,1} up to {Nx_max,Ny_max}. Using this
-- function, the table above can be simplified to:
--    return {
--      atmos = getmodes(2,2),
--      ocean = getmodes(2,2),
--    }
------------------------------------------------------------------------

--- Function to generate the mode blocks for either ocean or atmosphere up to
-- given wavenumbers.
-- @function getmodes
-- @param Nx_max Upper limit of the x wavenumber
-- @param Ny_max Upper limit of the y wavenumber
local function getmodes(Nx_max,Ny_max)
  local t = {}
  for Nx=1,Nx_max do
    for Ny=1,Ny_max do
      t[#t+1] = {Nx,Ny}
    end
  end
  return t
end

-- Atmospheric modes selection

-- Reinhold & Pierrehumbert 
local modes_RP82 = {{1,1},{1,2},{2,1},{2,2}} -- or: getmodes(2,2)
local modes_test = {{3,1},{3,2},{2,3},{1,1},{1,3},{1,2},{2,1},{2,2}
}

-- Ocean modes selection ( x block number accounts for half-integer wavenumber e.g 1 => 1/2 , 2 => 1, etc...) 

local modes_VD2014 = {{1,1},{1,2},{2,1},{2,2}}
local modes_VDDG14 = {{1,1},{1,2},{2,1},{2,2},{1,3},{1,4},{2,3},{2,4}}

return {atmos=modes_RP82,ocean=modes_VDDG14}
