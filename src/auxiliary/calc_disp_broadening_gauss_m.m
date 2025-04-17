function [broadening_rms,ld] = calc_disp_broadening_gauss_m(distance,beta2,params)
% Second-order dispersion-induced rms broadening of super-Gaussian pulses
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the root-mean-square broadening of a
% super-Gaussian pulse, possibly chirped, in a medium with second order
% dispersion described by its \beta_2 parameter.
% See for instance G. P. Agrawal, "Nonlinear fiber optics," chapter 3, 
% pp. 64–88, Academic Press, San Diego, 2nd edition, 1995.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params.fwhm = 10e-12;
% params.chirp = 0;
% params.order = 2;
% [broadening_rms,ld] = calc_disp_broadening_gauss_m(distance,beta2,params)
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% distance          distance at which the broadening will be calculated,
%                       in m [real vector]
%
% beta2             dispersion, in s^2/m [real scalar]
%
% params            super-Gaussian pulse parameters [structure]
%
%                       params.fwhm
%                           pulse full-width at half-maximum, in s
%                           [real scalar]
%
%                       params.chirp
%                           linear chirp parameter, no unit
%                           [real scalar]
%
%                       params.order
%                           Gaussian order [integer scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% broadening_rms    rms broadening factor [real vector]
%
% ld                dispersion length due to second-order dispersion, in m
%                       [real scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

t0 = params.fwhm/2/(log(2)^(1/2/params.order));
% Pulse t0 (half-width at 1/e)

ld = t0^2/beta2;
% Dispersion length, in m (signed)

distance = distance/ld;
% Normalisation of the propagation distance by the dispersion length

broadening_rms = sqrt(1 + (gamma(1/2/params.order)*params.chirp*distance)/gamma(3/2/params.order) +...
    (gamma(2 - 1/2/params.order)*(1 + params.chirp^2)*params.order^2*distance.^2)/gamma(3/2/params.order));
% Broadening factor
% See for instance G. P. Agrawal, "Nonlinear fiber optics," chapter 3, 
% pp. 64–88, Academic Press, San Diego, 2nd edition, 1995.

ld = abs(ld);
% Dispersion length, in m

end