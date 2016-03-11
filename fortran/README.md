# Modular arbitrary-order ocean-atmosphere model: MAOOAM -- Fortran implementation

------------------------------------------------------------------------

## About ##

(c) 2013-2016 Lesley De Cruz and Jonathan Demaeyer

See LICENSE.txt for license information.

This software is provided as supplementary material with:

* De Cruz, L., Demaeyer, J. and Vannitsem, S.: A modular arbitrary-order
ocean-atmosphere model: MAOOAM, Geosci. Model Dev. Discuss.,
[doi:xx/xxx](http://dx.doi.org/xx/xxx), 2016.

**Please cite this article if you use (a part of) this software for a
publication.**

The authors would appreciate it if you could also send a reprint of
your paper to <lesley.decruz@meteo.be>, <jonathan.demaeyer@meteo.be> and
<svn@meteo.be>. 

Consult the MAOOAM code repository at [doi:yy/yyy](http://dx.doi.org/yy/yyy)
for updates, and [our website](http://climdyn.meteo.be) for additional
resources.

------------------------------------------------------------------------

## Installation ##


The code should be compiled with gfortran 4.7+ (allows for allocatable arrays
in namelists). Unpack the archive in a folder, and run:
     make
 
 Remark: The command "make clean" removes the compiled files.

For Windows users, a minimalistic GNU development environment
 (including gfortran and make) is available at www.mingw.org .

------------------------------------------------------------------------

##  Description of the files ##

The model tendencies are represented through a tensor called aotensor which
includes all the coefficients. This tensor is computed once at the program
initialization.

* maooam.f90 : Main program.
* aotensor_def.f90 : Tensor aotensor computation module.
* IC_def.f90 : A module which loads the user specified initial condition.
* inprod_analytic.f90 : Inner products computation module.
* integrators.f90 : A module which contains the integrator for the model equations.
* Makefile : The Makefile.
* params.f90 : The model parameters module.
* README.txt : The present file.
* LICENSE.txt : The license text of the program.
* util.f90 : A module with various useful functions.
* stat.f90 : A module for statistic accumulation.
* params.nml : A namelist to specify the model parameters.
* int_params.nml : A namelist to specify the integration parameters.
* modeselection.nml : A namelist to specify which spectral decomposition will be used.
  
------------------------------------------------------------------------

## Usage ##

The user first has to fill the params.nml and int_params.nml namelist files according to their needs.

The modeselection.nml namelist can then be filled : 
* NBOC and NBATM specify the number of blocks that will be used in respectively the ocean and
  the atmosphere. Each block corresponds to a given x and y wavenumber.
* The OMS and AMS arrays are integer arrays which specify which wavenumbers of
  the spectral decomposition will be used in respectively the ocean and the
  atmosphere. Their shapes are OMS(NBOC,2) and AMS(NBATM,2).
* The first dimension specifies the number attributed by the user to the block and the second
  dimension specifies the x and the y wavenumbers.
* The VDDG model, described in Vannitsem et al. (2015) is given as an example
  in the archive.
* Note that the variables of the model are numbered according to the chosen
  order of the blocks.


Model and integration parameters can be specified in the params.nml namelist file.

Finally, the IC.nml file specifying the initial condition should be defined. To
obtain an example of this configuration file corresponding to the model you
have previously defined, simply delete the current IC.nml file (if it exists)
and run the program :

    ./maooam

It will generate a new one and start with the 0 initial condition. If you want another 
initial condition, stop the program, fill the newly generated file and restart :

    ./maooam

It will generate two files :
 * evol_field.dat : the recorded time evolution of the variables.
 * mean_field.dat : the mean field (the climatology)

------------------------------------------------------------------------

## Final Remarks ##

The authors would like to thank Kris for help with the lua2fortran project. It
has greatly reduced the amount of (error-prone) work.

  No animals were harmed during the coding process.
