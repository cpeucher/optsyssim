function hresp = calc_disp_imdd_small_signal(freq,lambda,dispersion,params)
% Fibre frequency response due to dispersion and chirp with direct detection
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the frequency response of a dispersive element
% due to cancellation of the beating terms between the carrier and
% symmetric sideband components upon direct detection in a photodiode. 
% The transfer function is calculated as a function of the modulation
% frequency and the chirp parameter (alpha) of the laser.
% The transfer function is considered as the ratio of the detected RF power
% after and before the dispersive element.
% The calculation is performed according to:
% J. Wang and K. Petermann, "Small signal analysis for dispersive optical 
% fiber communication systems," Journal of Lightwave Technology, vol. 10, 
% no. 1, pp. 96-100, January 1992 [doi: 10.1109/50.108743].
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% freq = [0:0.1:20]*1.0e9;                   % frequency range, in Hz
% lambda = CONSTANT.c/reference_frequency;   % laser wavelength, in m.
% dispersion = 17*80*1e-3;                 % dispersion, in s/m.
% params_chirpedtx.alpha = 0;                % laser alpha parameter.
% params_chirpedtx.fc = 0;                   % chirp corner frequency, in Hz.
% hresp = calc_disp_imdd_small_signal(freq,lambda,dispersion,params_chirpedtx); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% freq              frequency values at which the channel response will be
%                       evaluated, in Hz [real vector]
%
%
% lambda            laser wavelength, in m [real scalar]
%
% dispersion        fibre accumulated dispersion, in s/m
%                       [real scalar]
%
% params            laser chirp parameters [structure]
%
%                       params.alpha
%                           laser transient chirp (alpha) parameter
%
%                       params.fc
%                           laser adiabatic chirp corner frequency, in Hz
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% hresp             fibre transfer function (RF power after detection)
%                       [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% CONSTANT              essential physical constants [structure]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global CONSTANT

%--------------------------------------------------------------------------
% Calculation according to:
% F. Devaux, Y. Sorel, and J. F. Kerdiles, "Simple measurement of fiber 
% dispersion and of chirp parameter of  intensity modulated light emitter,"
% Journal of Lightwave Technology, vol. 11, no. 12, pp.1937-1940, Dec.
% 1993.
% Does not take adiabatic chirp into account.
%--------------------------------------------------------------------------

% hresp = sqrt(1 + params.alpha.^2)*cos(pi*lambda.^2*freq.^2*dispersion/CONSTANT.c + atan(params.alpha));
% % Fibre transfer function for light intensity / photocurrent.
% 
% hresp = hresp.*hresp;
% % The previous function is squared since the fiber transfer function is
% % defined here as the ratio of the RF power following photodetection (with
% % identical responsivity) after and before the dispersive element.

%--------------------------------------------------------------------------
% Calculation according to:
% Equation (33) in:
% J. Wang and K. Petermann, "Small signal analysis for dispersive optical 
% fiber communication systems," Journal of Lightwave Technology, vol. 10, 
% no. 1, pp. 96-100, January 1992 [doi: 10.1109/50.108743].
% Takes transient and adiabatic chirp into account
%--------------------------------------------------------------------------
theta = pi*lambda.^2*freq.^2*dispersion/CONSTANT.c;

hresp = abs(cos(theta)-sin(theta)*params.alpha.*(1 - 1i*params.fc./freq)).^2;
% Small signal frequency response.

end