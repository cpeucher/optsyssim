function osnr = meas_osnr(sig,params)
% Retrieve optical signal-to-noise ratio (OSNR) from an optical signal
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function retrieves the optical signal-to-noise ratio (OSNR) of an
% optical signal using the noise interpolation method. 
% It calculates the signal + noise power in a given measurement bandwidth
% around the signal frequency, then the noise power at a specified
% frequency offset on both sides of the signal. It then performs a linear
% interpolation of the noise spectral density at the signal frequency, from
% which it can retrieves the signal power from the signal + noise power
% calculation.
% This is a standard procedure for experimentally determining the OSNR
% using an optical spectrum analyser (OSA).
% Rectangular filters are used to select signal + noise at the signal 
% frequency as well as noise on both sides of the signal.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_osnr.sigfreq = params_tx.emission_frequency - reference_frequency;
% params_osnr.meas_bw = 100e9;
% params_osnr.meas_offset = 200e9;
% params_osnr.noise_bw = 12.5e9;
% osnr_retrieved_db = meas_osnr(sig,params_osnr);
%
% vline((params_osnr.sigfreq + params_osnr.meas_offset + params_osnr.meas_bw/2)/1.0e9,'r--')
% vline((params_osnr.sigfreq + params_osnr.meas_offset - params_osnr.meas_bw/2)/1.0e9,'r--')
% vline((params_osnr.sigfreq - params_osnr.meas_offset + params_osnr.meas_bw/2)/1.0e9,'r--')
% vline((params_osnr.sigfreq - params_osnr.meas_offset - params_osnr.meas_bw/2)/1.0e9,'r--')
% vline((params_osnr.sigfreq + params_osnr.meas_bw/2)/1.0e9,'b:')
% vline((params_osnr.sigfreq - params_osnr.meas_bw/2)/1.0e9,'b:')
% % Plot limits of measurement bandwidths on spectrum 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
% 
% params            noise interpolation measurement parameters [structure]
%
%                       params.sigfreq 
%                           signal centre frequency, relative to the 
%                           reference_frequency of the simulation, in Hz 
%                           [real scalar]
% 
%                       params.meas_bw
%                           bandwidth of the noise or signal + noise 
%                           measurement filter, in Hz [real scalar]
%
%                       params.meas_offset
%                           frequency offset from the signal at which the
%                           noise is measured, in Hz [real scalar]
%
%                       params.noise_bw      
%                           reference noise bandwidth in which the OSNR is
%                           calculated, in Hz [real scalar]
%
%                           Typically a value of 12.5 GHz (0.1 nm at 
%                           1550 nm) is used.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% osnr              retrieved optical signal-to-noise ratio, in dB
%                       [real scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% dt                time samples separation, in s [real scalar]
%
% frequency_array   relative frequency samples, in Hz [real vector]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global dt
global frequency_array

fmax = 1/dt/2;
% Maximum frequency of the optical frequency axis

if params.sigfreq + params.meas_offset + params.meas_bw/2 > fmax 
    error('meas_osnr: maximum frequency of upper sideband noise filter exceeds simulation bandwidth');
end

if params.sigfreq - params.meas_offset - params.meas_bw/2 < -fmax
    error('meas_osnr: minimum frequency of lower sideband noise filter exceeds simulation bandwidth');
end

if params.meas_bw > params.meas_offset 
    error('meas_osnr: signal and noise filters overlap');
end

params_obpf.type = 'rectangular_ideal';
params_obpf.centre_frequency = params.sigfreq + params.meas_offset;
params_obpf.bandwidth = params.meas_bw;
tf_p = opt_tf_obpf(params_obpf,frequency_array);
sig_p = opt_filter(sig,tf_p);
p_p = char_opt_average_power(sig_p);
% Upper sideband noise power measurement

params_obpf.centre_frequency = params.sigfreq - params.meas_offset;
tf_m = opt_tf_obpf(params_obpf,frequency_array);
sig_m = opt_filter(sig,tf_m);
p_m = char_opt_average_power(sig_m);
% Lower sideband noise power measurement

params_obpf.centre_frequency = params.sigfreq;
tf_0 = opt_tf_obpf(params_obpf,frequency_array);
sig_0 = opt_filter(sig,tf_0);
p_0 = char_opt_average_power(sig_0);
% Signal + noise power measurement

osnr = 10*log10((p_0 - 0.5*(p_m + p_p))/(0.5*(p_m + p_p)*params.noise_bw/params_obpf.bandwidth));
% OSNR calculation by linear interpolation of the noise PSD at the signal
% frequency

end