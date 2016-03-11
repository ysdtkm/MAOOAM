-- gz.lua
-- (C) 2015-2016 Lesley De Cruz
-- See LICENSE.txt for license information.

------------------------------------------------------------------------
-- gz.lua is a very basic ffi binding to the functions gzopen, gzwrite and
-- gzclose of zlib.
------------------------------------------------------------------------

local ffi = require("ffi")
ffi.cdef
[[
typedef void *voidp;
typedef void const *voidpc;
typedef voidp gzFile;
typedef struct { gzFile f;} gzf;

gzFile
gzopen (const char *path , const char *mode );
int
gzwrite (gzFile file, voidpc buf, unsigned int len);
int
gzclose (gzFile file);
]]

local zlib = ffi.load(ffi.os == "Windows" and "zlib1" or "z")

--- Wrapper around gzopen
-- @param fname string containing filename
-- @param mode string defining the open mode (read: "r", write: "w", binary:
-- "b") Default: "wb"
-- @return a gzf object with write and close methods.
local function gzopen(fname, mode)
  return ffi.new("gzf",zlib.gzopen(fname,mode or "wb"))
end

--- Wrapper around gzwrite
-- @param gzf a gzf object (the result of gzopen)
-- @param str a string to compress and write to gzf
-- @return the number of processed bytes (useful for assertions)
local function gzwrite(gzf,str)
  str = tostring(str)
  return zlib.gzwrite(gzf.f,str,#str)
end

--- Wrapper around gzclose.
-- @param gzf a gzf object (the result of gzopen)
local function gzclose(gzf)
  return zlib.gzclose(gzf.f)
end

-- @export
return ffi.metatype("gzf",{__index={open = gzopen, write = gzwrite, close = gzclose}})
