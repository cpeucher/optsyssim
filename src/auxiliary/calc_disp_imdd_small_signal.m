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
% The transfer function is considered as the ratio of the detected current
% after and before the dispersive element.
% It corresponds to the S21 parameter one would measure with a vector
% network analyser.
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
% dispersion = 17*80*1e-3;                   % dispersion, in s/m.
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

theta = pi*lambda.^2*freq.^2*dispersion/CONSTANT.c;
% theta = -0.5*beta_2 w^2 L, where w = angular frequency

hresp = cos(theta)-sin(theta)*params.alpha.*(1 - 1i*params.fc./freq);
% Small signal frequency response
% Ratio of currents / voltage => S21
% In dB: 20*log10(abs(hresp))

end