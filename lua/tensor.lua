-- tensor.lua
-- (C) 2013 Lesley De Cruz
-- See LICENSE.txt for license information.

------------------------------------------------------------------------
-- Sparse tensor tools.
--
-- Tensors are entered as coordinate lists:
--    { {x1,y1,z1,val1}, {x2,y2,z2,val2},...}
--
-- This module provides functions to convert a coordinate lists to a nested
-- sparse table or to a (non-growable) fast ffi coordinate list.
--
-- While not the most (space-)efficient sparse tensor storage method, it's
-- simple, generic (for variable dimensions) and fast enough for
-- contractions with 1D arrays.
------------------------------------------------------------------------

local ffi = require("ffi")
local ipairs, tonumber, assert, yield = ipairs, tonumber, assert, coroutine.yield
if not DEBUG then assert = function(v) return v end end

--- Convert a coordinate list of arbitrary dimension to a sparse nested table.
-- Input format: { {x1,y1,z1,val1}, {x2,y2,z2,val2},...}
-- In the returned table, we have
--     t[x1][y1][z1]==val1
-- e.g.
--     mycoolist = {{1, 5, 0.1},{2, 3, 3.14}}
-- yields
--     {[1]={[5]=0.1},[2]={[3]=3.14}}
-- Missing entries are treated as zero in mathematical operations.
-- @function coo_to_sparse_table
-- @param coo_list coordinate list (Lua table)
-- @return sparse table
local coo_to_sparse_table

local new_sparse_table

do
  local sparse_mt
  local function asnumber(x) return tonumber(x) or 0 end
  sparse_mt = {
    -- Assignment function for sparse tables.
    assign = function(sparse_t,value,...)
      local narg = select("#",...)
      for i=1,narg-1 do
        local j = select(i,...)
        sparse_t[j]= rawget(sparse_t,j) or setmetatable({},sparse_mt)
        sparse_t = sparse_t[j]
      end
      sparse_t[select(narg,...)] = value
    end,
    -- Automagically generate tables. 
    __index = function(t,k) 
      if tonumber(k) then return setmetatable({},sparse_mt)
      else return sparse_mt[k] end 
    end,
    -- Nonexistent entries behave as zero in mathematical operations.
    __add = function(a,b) return asnumber(a) + asnumber(b) end,
    __sub = function(a,b) return asnumber(a) - asnumber(b) end,
    __unm = function(a) return 0 end,
    __mul = function(a,b) return asnumber(a) * asnumber(b) end,
    __div = function(a,b) return asnumber(a) / asnumber(b) end,
  }
  -- Recursively insert a value dim coordinates deep into a table.
  local function insert_nested(t,level,dim,coo)
    if level==dim then
      -- If more than one entry at a given set of coordinates, take the sum.
      t[coo[level]] = t[coo[level]] + coo[dim+1]
      return
    end
    t[coo[level]] = t[coo[level]] or setmetatable({},sparse_mt)
    insert_nested(t[coo[level]],level+1,dim,coo)
  end
  function coo_to_sparse_table(coo_list)
    local t = setmetatable({},sparse_mt)
    local dim = #coo_list[1]-1
    for _,coo in ipairs(coo_list) do
      insert_nested(t,1,dim,coo)
    end
    return t
  end
  --- Create a new sparse table.
  -- Missing entries are treated as zero in mathematical operations.
  -- This sparse table cannot be filled up with st[i][j][k]=val
  -- Instead, use the assign method sparse_t:assign(val,indices...)
  -- @function new_sparse_table
  -- @return sparse table
  function new_sparse_table() return setmetatable({},sparse_mt) end
end

-- Get the next entry from a nested sparse table in coordinate list form:
-- val1, z1, y1, x1
-- Note that the entries are reversed!
local function next_coo(t,...)
  if type(t)~="table" then yield(t,...)
  else
    for k,t2 in pairs(t) do next_coo(t2,k,...) end
  end
end

-- Reverse the entries of a table.
local function treverse(t)
  for k=1,#t/2 do
    t[k],t[#t+1-k]=t[#t+1-k],t[k]
  end
  return t
end

--- Convert a sparse tensor encoded as a nested table into a coordinate list.
-- Don't include elements with value equal to zero.
-- @param sparse_table sparse, nested table
-- @return coordinate list (Lua table)
local function sparse_table_to_coo(sparse_table)
  local getnextcoo = coroutine.wrap(function() return next_coo(sparse_table) end)
  local coo_list = {}
  while true do
    local coo = treverse{getnextcoo()}
    if coo[1] then
      if coo[#coo]~=0 then
        coo_list[#coo_list+1] = coo
      end
    else break end
  end
  return coo_list
end

--- Merge duplicate entries in a coordinate list.
-- This is done by converting it to a sparse table and back.
-- @param coo_list coordinate list (Lua table)
-- @return simplified coordinate list (Lua table)
local function simplify_coo(coo_list)
  return sparse_table_to_coo(coo_to_sparse_table(coo_list))
end

local fficoo_cache = {}
--- Convert a coordinate list of arbitrary dimension to an ffi coordinate list.
-- Input format:
--     { {x1,y1,z1,val1}, {x2,y2,z2,val2},...}
-- e.g.
--     mycoolist = {{1, 5, 0.1}, {2, 3, 3.14}, {3, 2, 10}}
-- yields a struct fficoo which has
--     fficoo.dim == 2
--     fficoo.nelem == 3
--     fficoo.data[0].c[0] == 1
--     fficoo.data[0].c[1] == 5
--     fficoo.data[0].v == 0.1
--     fficoo.data[1].c[0] == 2
--  etc.
-- @param coo_list coordinate list (Lua table)
-- @return ffi coordinate list (ffi struct)
local function coo_to_fficoo(coo_list)
  local nelem = #coo_list
  assert(nelem>0, "coo_to_fficoo can't handle an empty coo_list!")
  local dim = #coo_list[1]-1
  fficoo_cache[dim] = fficoo_cache[dim] or {}
  local fficoo_t = fficoo_cache[dim][nelem]
  if not fficoo_t then
    local coo_entry_t = ffi.typeof("struct {double v; int32_t c[$];}",dim)
    fficoo_t = ffi.typeof("struct{ int32_t dim,nelem; $ data[$];}",coo_entry_t,nelem)
    fficoo_cache[dim][nelem] = fficoo_t
  end
  local inittab = {dim=dim, nelem=nelem, data={}} -- table used to initialize the fficoo struct.
  for i,coo in ipairs(coo_list) do
    local init_coo = {v=coo[dim+1], c={}}
    for j=1,dim do init_coo.c[j-1]=coo[j] end
    inittab.data[i] = init_coo
  end
  return fficoo_t(inittab)
end

--- Sparse multiplication of two tensors, C[i][j] a[j].
-- Note that it is NOT safe to pass `arr_j` as a result buffer, 
-- as this operation does multiple passes.
-- @param fficoo_ij an ffi coordinate list (ffi struct) of which index
-- 2 will be contracted.
-- @param arr_j the array to be contracted with index 2 of ffi_coo_ij
-- @param res array (buffer) to store the result of the contraction
-- @return the result array
local function sparse_mul2(fficoo_ij, arr_j, res)
  res:clear()
  local dim = fficoo_ij.dim
  assert(dim==2, "sparse_mul2 expects a coordinate list with exactly 2 coordinates")
  assert(fficoo_ij.data[1].v, "sparse_mul2 expects an ffi coordinate list.")
  for ci = 0,fficoo_ij.nelem-1 do -- sum over j
    local coo = fficoo_ij.data[ci]
    res.data[coo.c[0]] = res.data[coo.c[0]] + coo.v*arr_j.data[coo.c[1]]
  end
  return res
end

--- Sparse multiplication of three tensors, C[i][j][k] a[j] b[k].
-- Note that it is NOT safe to pass `arr_j`/`arr_k` as a result buffer, 
-- as this operation does multiple passes.
-- @param fficoo_ijk an ffi coordinate list (ffi struct) of which index
-- 2 and 3 will be contracted.
-- @param arr_j the array to be contracted with index 2 of ffi_coo_ijk
-- @param arr_k the array to be contracted with index 3 of ffi_coo_ijk
-- @param res array (buffer) to store the result of the contraction
-- @return the result array
local function sparse_mul3(fficoo_ijk, arr_j, arr_k, res)
  res:clear()
  local dim = fficoo_ijk.dim
  assert(dim==3, "sparse_mul3 expects a coordinate list with exactly 3 coordinates")
  assert(fficoo_ijk.data[1].v, "sparse_mul3 expects an ffi coordinate list.")
  for ci = 0,fficoo_ijk.nelem-1 do -- sum over j and k
    local coo = fficoo_ijk.data[ci]
    res.data[coo.c[0]] = res.data[coo.c[0]] + coo.v*arr_j.data[coo.c[1]]*arr_k.data[coo.c[2]]
  end
  return res
end

--- Sparse multiplication of two tensors to determine the Jacobian:
--    J[i][j] = C[i][j][k] a[k] + C[i][k][j] a[k].
-- It's implemented slightly differently: for every C[i][j][k], we add to J as follows:
--    J[i][j] += C[i][j][k] a[k];
--    J[i][k] += C[i][j][k] a[j]
--- Return a tensor in the form of a coolist.
---- @param fficoo_ijk an ffi coordinate list (ffi struct) of which index
---- 2 or 3 will be contracted.
---- @param arr_j the array to be contracted with index 2 and then index 3 of ffi_coo_ijk
---- @return jcoo_ij a coolist with the result of the contraction
local function jsparse_mul(fficoo_ijk, arr_j)
  local jcoo_ij={}
  local dim = fficoo_ijk.dim
  assert(dim==3, "jsparse_mul expects a coordinate list with exactly 3 coordinates")
  assert(fficoo_ijk.data[1].v, "jsparse_mul expects an ffi coordinate list.")
  for ci = 0,fficoo_ijk.nelem-1 do -- sum over j and k
    local coo = fficoo_ijk.data[ci]
    if coo.c[1] ~= 0 then jcoo_ij[#jcoo_ij+1]= {coo.c[0],coo.c[1],coo.v*arr_j.data[coo.c[2]]} end
    if coo.c[2] ~= 0 then jcoo_ij[#jcoo_ij+1]= {coo.c[0],coo.c[2],coo.v*arr_j.data[coo.c[1]]} end
  end
  return jcoo_ij
end


--- Sparse multiplication of four tensors, C[i][j][k][l] a[j] b[k] c[l]
-- Note that it is NOT safe to pass `arr_j/k/l` as a result buffer, 
-- as this operation does multiple passes.
-- @param fficoo_ijkl an ffi coordinate list (ffi struct) of which index
-- 2,3,4 will be contracted.
-- @param arr_j the array to be contracted with index 2 of ffi_coo_ijkl
-- @param arr_k the array to be contracted with index 3 of ffi_coo_ijkl
-- @param arr_l the array to be contracted with index 4 of ffi_coo_ijkl
-- @param res array (buffer) to store the result of the contraction
-- @return the result array
local function sparse_mul4(fficoo_ijkl, arr_j, arr_k, arr_l, res)
  res:clear()
  local dim = fficoo_ijkl.dim
  assert(dim==4, "sparse_mul4 expects a coordinate list with exactly 4 coordinates")
  assert(fficoo_ijkl.data[1].v, "sparse_mul4 expects an ffi coordinate list.")
  for ci = 0,fficoo_ijkl.nelem-1 do -- sum over j and k
    local coo = fficoo_ijkl.data[ci]
    res.data[coo.c[0]] = res.data[coo.c[0]] + coo.v*
     arr_j.data[coo.c[1]]*arr_k.data[coo.c[2]]*arr_l.data[coo.c[3]]
  end
  return res
end

--- Sparse multiplication of five tensors, C[i][j][k][l][m] a[j] b[k] c[l] d[m]
-- Note that it is NOT safe to pass `arr_j/k/l/m` as a result buffer, 
-- as this operation does multiple passes.
-- @param fficoo_ijklm an ffi coordinate list (ffi struct) of which index
-- 2,3,4,5 will be contracted.
-- @param arr_j the array to be contracted with index 2 of ffi_coo_ijklm
-- @param arr_k the array to be contracted with index 3 of ffi_coo_ijklm
-- @param arr_l the array to be contracted with index 4 of ffi_coo_ijklm
-- @param arr_m the array to be contracted with index 5 of ffi_coo_ijklm
-- @param res array (buffer) to store the result of the contraction
-- @return the result array
local function sparse_mul5(fficoo_ijklm, arr_j, arr_k, arr_l, arr_m, res)
  res:clear()
  local dim = fficoo_ijklm.dim
  assert(dim==5, "sparse_mul5 expects a coordinate list with exactly 5 coordinates")
  assert(fficoo_ijklm.data[1].v, "sparse_mul5 expects an ffi coordinate list.")
  for ci = 0,fficoo_ijklm.nelem-1 do -- sum over j and k
    local coo = fficoo_ijklm.data[ci]
    res.data[coo.c[0]] = res.data[coo.c[0]] + coo.v*
     arr_j.data[coo.c[1]]*arr_k.data[coo.c[2]]*arr_l.data[coo.c[3]]*arr_m.data[coo.c[4]]
  end
  return res
end

--- Append a coo_list to another one.
-- @param coo1 the target coo_list
-- @param coo2 the coo_list that will be appended to coo1
local function append_coo(coo1,coo2)
  for k,v in ipairs(coo2) do
    coo1[#coo1+1]=v
  end
end

--- @export
return {
  coo_to_sparse_table = coo_to_sparse_table,
  new_sparse_table = new_sparse_table,
  sparse_table_to_coo = sparse_table_to_coo,
  simplify_coo = simplify_coo,
  coo_to_fficoo = coo_to_fficoo,
  sparse_mul2 = sparse_mul2,
  sparse_mul3 = sparse_mul3,
  sparse_mul4 = sparse_mul4,
  sparse_mul5 = sparse_mul5,
  jsparse_mul = jsparse_mul,
  append_coo = append_coo,
}
