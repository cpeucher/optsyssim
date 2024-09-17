function sig = tx_duobinary_nrz_low_pass(params)
% Low-pass filter NRZ duobinary transmitter
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates an ideal optical duobinary signal by filtering 
% the 2-level electrical NRZ signal with a Bessel low pass filter with 
% bandwidth of 1/4 of the bit rate, before applying it to a Mach-Zehnder 
% modulator biased at a transmission minimum.
% The transmitter consists of a CW laser followed by a Mach-Zehnder data 
% modulator. The generated signal is ideal in the sense that the extinction
% ratios of the data modulator is infinite and the CW laser has zero 
% linewidth.
% The data modulator is biased at minimum point and driven with a 
% peak-to-peak voltage difference between the 2 driving signals equal to 
% 2Vpi.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_tx.emission_frequency = reference_frequency;
% params_tx.power = 1.0e-3;
% params_tx.linewidth = 0;
% params_tx.bit_pattern = bit_pattern;
% params_tx.rise_time = 1/symbol_rate/4;
% sig = tx_duobinary_nrz_low_pass(params_tx); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% params            transmitter parameters [structure]
%
%                       params.emission_frequency
%                           desired centre frequency of the channel, in Hz
%                           [real scalar]
%                           It will be adapted so that it coincides with a
%                           sample of the discrete frequency grid.
%
%                       params.linewidth
%                           laser linewidth, in Hz [real scalar]
%
%                       params.laser_power
%                           laser power, in W [real scalar]
%
%                       params.bit_pattern
%                           binary data pattern of the signal to 
%                           be generated [binary vector]
%
%                       params.rise_time
%                           rise time of the electrical NRZ signal driving
%                           the data modulator, in s [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               modulated optical signal [optical signal structure]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% dt                    time samples separation, in s [real scalar]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global dt


bit_pattern_diff_enc = logical_differential_encoder_binary(not(params.bit_pattern));
% Pre-coding

nrz_data_sig = elec_pulse_sequence_nrz(bit_pattern_diff_enc,params.rise_time);
% Generate differentially-encoded NRZ data stream


bit_rate = length(params.bit_pattern)/length(nrz_data_sig)/dt;
% Signal bit rate

elpf.type = 'bessel';
elpf.order = 4;
elpf.f3dB = 0.25*bit_rate;
duobinary_data_sig = elec_elpf(nrz_data_sig,elpf);
% Low pass filtering for duobinary encoding

sig = opt_laser_cw(params);
% CW laser signal

vpi = 1.0;
driving_signal_1 = vpi*(duobinary_data_sig - 0.5);
driving_signal_2 = -vpi*(duobinary_data_sig - 0.5);
bias_1 = vpi;
bias_2 = 0;
sig = mod_mzm(sig,driving_signal_1,driving_signal_2,bias_1,bias_2,vpi,0.5,0.5,0);
% Data modulator

[sig,~] = rx_resynchronise(sig,params.bit_pattern,'opt');
% Resynchronise output signal to compensate for possible delay

end