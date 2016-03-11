-- rk2.lua
-- (C) 2013 Lesley De Cruz
-- See LICENSE.txt for license information.

------------------------------------------------------------------------
-- A simple yet stable second order Runge-Kutta scheme (Heun's method).
-- Butcher tableau:
--
--    0 |
--    1 | 1
--      +----------
--        1/2  1/2
------------------------------------------------------------------------

--- Create an integrator for n-arrays.
-- @function rk2
-- @param n number of variables (length of the array)
-- @return Integrator for n-arrays
return function(n)
  local array = require("array")(n)
  local buf_f0, buf_f1, buf_y1 = array(), array(), array()
  --- Integrator for n-arrays
  -- @function integrator
  -- @param y variables at time t
  -- @param f function to calculate the time derivatives of the variables
  -- @param t time
  -- @param dt time integration step
  -- @param ynew n-array (buffer) to store the new y-value
  -- @return t+dt (incremented time)
  -- @return ynew array with variables at time t+dt
  return function(y,f,t,dt,ynew)
    f(t,y,buf_f0)                        -- f0 = f(t,y)
    buf_f0:nmul(dt,buf_y1):add(y,buf_y1)  -- y1 = y + f0*dt
    f(t+dt,buf_y1,buf_f1)                -- f1 = f(t+dt, y1)
    -- ynew = y + (f0+f1)*0.5*dt
    return t+dt, buf_f0:add(buf_f1,buf_f0):nmul(0.5*dt,buf_f0):add(y,ynew)
  end
end

