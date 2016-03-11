-- rk4.lua
-- (C) 2015 Lesley De Cruz
-- See LICENSE.txt for license information.

------------------------------------------------------------------------
-- A classical fourth order Runge-Kutta scheme. 
-- Butcher tableau:
--
--    0   |
--    1/2 | 1/2
--    1/2 | 0    1/2
--    1   | 0    0    1
--        +-------------------
--          1/6  1/3  1/3  1/6
------------------------------------------------------------------------

--- Create an integrator for n-arrays.
-- @function rk4
-- @param n number of variables (length of the array)
-- @return Integrator for n-arrays
return function(n)
  local array = require("array")(n)
  local buf_kA, buf_kB, buf_y1 = array(), array(), array()
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
    f(t,y,buf_kA)                            -- k_1 (A) = f(t, y)
    buf_kA:nmul(0.5*dt,buf_y1):add(y,buf_y1) -- y_1     = y + 0.5*dt*k_1 (A)
    f(t+0.5*dt,buf_y1,buf_kB)                -- k_2 (B) = f(t + 0.5*dt, y_1)
    buf_kB:nmul(0.5*dt,buf_y1):add(y,buf_y1) -- y_1'    = y + 0.5*dt*k_2
    buf_kB:nmul(2,buf_kB):add(buf_kA,buf_kA) -- k_S (A) = k_1 (A) + 2 * k_2 (B)
    f(t+0.5*dt,buf_y1,buf_kB)                -- k_3 (B) = f(t + 0.5*dt, y_1')
    buf_kB:nmul(dt,buf_y1):add(y,buf_y1)     -- y_1''   = y + dt*k_3 (B)
    buf_kB:nmul(2,buf_kB):add(buf_kA,buf_kA) -- k_S (A) = k_S (A) + 2 * k_3 (B)
    f(t+dt,buf_y1,buf_kB)                    -- k_4 (B) = f(t + dt, y_1")
    buf_kB:add(buf_kA,buf_kA)                -- k_S (A) = k_S (A) + k_4 (B)
     -- ynew = y + dt/6 * k_S
    return t+dt, buf_kA:nmul(dt/6,buf_kA):add(y,ynew)
  end
end

