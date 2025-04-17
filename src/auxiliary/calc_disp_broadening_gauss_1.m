function [broadening_rms,ld,ldp] = calc_disp_broadening_gauss_1(distance,beta2,beta3,params)
% 2nd and 3rd order dispersion-induced rms broadening of 1st order Gaussian pulse
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the root-mean-square broadening of a
% 1st order Gaussian pulse, possibly chirped, under the influence of second
% and third order dispersion.
% See for instance G. P. Agrawal, "Nonlinear fiber optics," chapter 3, 
% pp. 64–88, Academic Press, San Diego, 2nd edition, 1995.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params.fwhm = 10e-12;
% params.chirp = 0;
% [broadening_rms,ld,ldp] =
% calc_disp_broadening_gauss_1(distance,beta2,beta3,params);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% distance          distances at which the broadening will be calculated,
%                       in m [real vector]
%
% beta2             second-order dispersion, in s^2/m [real scalar]
%
% beta3             third-order dispersion, in s^3/m [real scalar]
%
% params            1st order Gaussian pulse parameters [structure]
%
%                       params.fwhm
%                           pulse full-width at half-maximum, in s
%                           [real scalar]
%
%                       params.chirp
%                           linear chirp parameter, no unit
%                           [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% broadening_rms    rms broadening factor, no unit [real vector]
%
% ld                dispersion length due to second-order dispersion, in m
%                       [real scalar]
%
% ldp               dispersion length due to third-order dispersion, in m
%                       [real scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

t0 = params.fwhm/2/(log(2)^(1/2));
% Pulse t0 (half-width at 1/e)
% Linked to the rms duration according to 
% rms0 = t0/sqrt(2)

ld = t0^2/beta2;
% Dispersion length due to 2nd order dispersion, in m (signed)
ldp = t0^3/beta3;
% Dispersion length due to 3rd order dispersion, in m (signed)

broadening_rms = sqrt((1 + params.chirp*distance/ld).^2 + (distance/ld).^2 + 0.25*(1 + params.chirp^2)^2*(distance/ldp).^2);
% Broadening factor
% See for instance G. P. Agrawal, "Nonlinear fiber optics," chapter 3, 
% pp. 64–88, Academic Press, San Diego, 2nd edition, 1995.

ld = abs(ld);
% Dispersion length linked to 2nd order dispersion, in m

ldp = abs(ldp);
% Dispersion length linked to 3rd order dispersion, in m

end