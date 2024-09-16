function [sig,freq_data] = tx_laser_chirped(sig,params)
% Black-box laser with frequency chirp
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements a simplified black-box model of a laser
% described by the characteristics of its emitted signal, namely average
% power, extinction ratio, alpha parameter and adiabatic chirp parameter.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_chirped_tx.emission_frequency = reference_frequency;   % laser emission frequency, in Hz
% params_chirped_tx.power_average = 1.0e-3;                     % signal average power, in W
% params_chirped_tx.extinction_ratio = 10;                      % signal extinction ratio, in dB
% params_chirped_tx.alpha = 2.5;                                % alpha parameter
% params_chirped_tx.kappa = 12e12;                              % adiabatic chirp, in Hz/W
% [sig,freq_data] = tx_laser_chirped(nrz_data_sig,params_chirped_tx); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input electrical signal [real vector]
%
% params            signal parameters [structure]
%
%                       params.emission_frequency
%                           laser emission frequency, in Hz [real scalar]
%
%                       params.power_average
%                           signal average power, in W [real scalar]
%
%                       params.extinction_ratio
%                           signal extinction ratio, in dB [real scalar]
%
%                       params.alpha
%                           transcient chirp alpha parameter [real scalar]
%
%                       params.kappa
%                           diabatic chirp parameter, in Hz/W [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               modulated optical signal [optical signal structure]
%
% freq_data         information about the emission frequencies for the 
%                       marks and spaces [structure]
%
%                       freq_data.nu0
%                           emission frequency of the spaces, in Hz
%                           [real scalar]
%
%                       freq_data.nu1
%                           emission frequency of the marks, in Hz
%                           [real scalar]
%
%                       freq_data.frequency_offset
%                           frequency offset due to phase discontinuity 
%                           compensation, in Hz [real scalar]
%
%                       freq_data.nu0 and freq_data.nu1 are expressed 
%                       relatively to the emission frequency of the
%                       laser described by the input parameter 
%                       params.emission_frequency and do not take the 
%                       frequency offset due to phase discontinuity 
%                       compensation into account.
%                       They are are result of analytical calculations
%                       of the adiabatic chirp, assuming that the
%                       driving signal takes values within [0 1], and
%                       reaches the limits of this interval, i.e. is
%                       not excessively low-pass filtered.
%                       The actual emission frequencies for the spaces
%                       and marks are:
%                       params.emission_frequency + freq_data.nu0 + freq_data.frequency_offset
%                       and
%                       params.emission_frequency + freq_data.nu0 + freq_data.frequency_offset
%                       respectively.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% time_array            time samples, in s [real vector]
%
% dt                    time samples separation, in s [real scalar]
%
% reference_frequency   reference frequency, in Hz [real scalar]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global dt 
global reference_frequency 
global time_array


sig = sig(:).';
% Ensure the input signal is a line vector.

nsamples = length(sig);
% Number of samples in the electrical input signal.

sig = (sig - min(sig)) / (max(sig) - min(sig));
% Ensure the modulating signal takes values within [0 1].

extinction_ratio_linear = 10^(params.extinction_ratio/10);
% Convert extinction ratio from logarithmic to linear scale.

power_peak_peak = 2*params.power_average*(extinction_ratio_linear - 1)/(extinction_ratio_linear + 1);
% Peak-to-peak power of the emitted signal.

signal_power = params.power_average - power_peak_peak/2+power_peak_peak*sig;
% Power as  a function of time.

signal_phase = -0.5*params.alpha*((log(signal_power/signal_power(1)) + params.kappa*cumtrapz(signal_power)*dt));
% Phase as a function of time, taking adiabatic and transient chirp into
% account.


% Analytical calculation of the frequencies of the marks and spaces,
% assuming that the driving signal is not excessively low-pass filtered.
% These values will be returned, together with the frequency offset due to
% compensation of phase discontinuities (see below) in order to assist the
% precise frequency tuning of the laser.
power_spaces = 2*params.power_average/(extinction_ratio_linear + 1);
% Power of the spaces.
power_marks = extinction_ratio_linear*power_spaces;
% Power of the marks.
freq_data.nu0 = params.alpha*params.kappa*power_spaces/4/pi;
% Frequency of the spaces.
freq_data.nu1 = params.alpha*params.kappa*power_marks/4/pi;
% Frequency of the marks.

% Apply correction factor to the phase noise process in order to avoid
% phase discontinuities in conjunction with the periodic nature of the
% signal.
% C. J. Rasmussen, "Transmission analysis in WDM networks", Ph.D. thesis,
% Technical University of Denmark, 1999.
% pp. 118-119.
phase_diff = signal_phase(nsamples) - signal_phase(1);
if phase_diff ~= 0
    % No need to worry if the linewidth is zero or if the first and last
    % samples have the same phase (the phase may be continuous though with
    % PhaseDiff ~= 0 since phase(1)=phase(Nsamples +1) in this case.     
    delta_phi = round(phase_diff/(2*pi))*2*pi;
    phase_correction = time_array*(delta_phi-phase_diff)/((nsamples - 1)*dt);
    % Phase correction to ensure continuity modulo 2*pi of the phase.
    
    freq_data.frequency_offset = -(delta_phi - phase_diff)/((nsamples - 1)*dt)/2/pi;
    % Frequency offset due to phase discontinuity compensation.   
    
    signal_phase = signal_phase + phase_correction;
    % Apply phase correction.
else
    freq_data.frequency_offset = 0;
    % Frequency offset due to phase discontinuity compensation.
    
end

sig = struct;
% Initialise optical signal structure.

sig.x = sqrt(signal_power).*exp(-1i*signal_phase).*exp(1i*2*pi*(params.emission_frequency - reference_frequency)*time_array);
sig.y = zeros(1,length(sig));
% Field (complex envelope with respect to reference_frequency) of the
% optical emitted signal, calculated from the knowledge of its power, phase
% and emission frequency.


% power_average_retrieved = char_opt_average_power(sig);
% Calculate average power from the signal itself.

% extinction_ratio_estimated = 10*log10(max(abs(sig.x).^2)/min(abs(sig.x).^2));
% Estimate of the extinction ratio of the signal. This is valid only for
% "clean" driving signals that have not bee excessively low-pass filtered.


end