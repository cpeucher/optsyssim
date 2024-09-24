function CONSTANT = core_load_constants()
% Load essential physical constants
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function loads essential physical constants of interest. 
% The values are based on CODATA Recommended values 2022.
% see e.g. http://physics.nist.gov/cuu/Constants/index.html
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% CONSTANT = core_load_constants();
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% None
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% CONSTANT              essential physical constants [structure]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

CONSTANT.q = 1.602176634e-19;
% Elementary charge, in C
CONSTANT.c = 299792458;
% Speed of light in vacuum, in m/s
CONSTANT.h = 6.62607015e-34;
% Planck constant, in J.s
CONSTANT.hbar = CONSTANT.h/2/pi;
% Planck's constant (hbar), in J.s
CONSTANT.kb = 1.380649e-23;
% Boltzmann's constant, in J/K
CONSTANT.eps0 = 8.8541878188e-12;
% Vacuum electric permittivity, in F/m
CONSTANT.mu0 = 1.25663706127e-6;
% Vacuum magnetic permeability, in N/A^2
CONSTANT.Z0 = 376.730313412;
% Characteristic impedance of vacuum, in ohm
CONSTANT.me = 9.1093837139e-31;
% Electron mass, in kg




end