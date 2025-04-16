function [rf_power,psd] = meas_rf_power(sig,params)
% Integrate RF power in a specified bandwidth
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function determines the RF power of an electrical signal in a given 
% bandwidth around a given centre frequency.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_rf_pwm.centre_frequency = 10e9;
% params_rf_pwm.bandwidth = 4*df;
% params_rf_pwm.input_impedance = 1;
% [rf_power,psd] = meas_rf_power(sig,params_rf_pwm);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               electrical signal to be characterised [real vector]
%                       
% params            structure containing RF power measurement parameters
%                       [structure]
%
%                       params.centre_frequency
%                           centre frequency around which the RF power will 
%                           be integrated [real scalar]
%
%                       params.bandwidth
%                           rectangular bandwidth over which the RF power 
%                           will be integrated [real scalar]
%                            Better to be an integer multiple of df.
%
%                       params.input_impedance
%                           input resistance for calculation of the RF 
%                           power, in ohms [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% rf_power          RF power in params.bandwidth around 
%                       params.centre_frequency, in W [real scalar]
%
% psd               estimated power spectral density, in W/Hz [real scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% frequency_array   relative frequency samples, in Hz [real vector]
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global frequency_array


% -------------------------------------------------------------------------
% Check frequency range over which power will be integrated
% -------------------------------------------------------------------------
if params.centre_frequency - params.bandwidth/2 < 0
    error('meas_rf_power: desired bandwidth too large for the desired centre frequency.');
end

% -------------------------------------------------------------------------
% Calculate one-sided spectrum
% -------------------------------------------------------------------------
nsamples = length(sig);
% Number of samples in the signal

xf = fft(sig)/nsamples;
% fft of the input signal, in fft order

xf = xf(1:nsamples/2)/sqrt(params.input_impedance);
% Restrict to positive frequency range (incl. dc) and take load resistance
% into account
xf(2:end) = sqrt(2)*xf(2:end);
% Multiply by sqrt(2) apart from the DC term, so that the single-sided 
% spectrum is considered
% xf contains the one-sided spectrum over the positive frequencies,
% including dc

freq_rf = frequency_array(nsamples/2+1:end);
% Positive frequencies


% -------------------------------------------------------------------------
% Filtering over specified bandwidth and integration of the power
% -------------------------------------------------------------------------
tf = abs(freq_rf - params.centre_frequency) <= params.bandwidth/2;   
% Transfer function of the rectangular filter used for bandwidth limitation
% to the specified bandwidth around params.centre_frequency

xf = xf.*tf;
% Applies transfer function to the signal spectrum.
% OBS: this is the Fourier transform of the signal. We do not calculate the
% power spectrum here

rf_power = sum(abs(xf).^2);
% Total power integrated in the specified bandwidth around
% centre_frequency

psd = rf_power/params.bandwidth;
% Corresponding power spectral density

end