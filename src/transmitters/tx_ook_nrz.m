function sig = tx_ook_nrz(params)
% NRZ-OOK transmitter
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates an ideal NRZ-OOK signal.
% The transmitter consists of a CW laser followed by a Mach-Zehnder data 
% modulator. The generated signal is ideal in the sense that the extinction
% ratios of the data modulator is infinite, and chirp-free modulation is 
% applied.
% The data modulator is biased at quadrature point and driven with a 
% peak-to-peak voltage difference between the 2 driving signals equal to 
% Vpi.
% For other signal generation configurations, directly use the mod_mzm
% model.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_tx.emission_frequency = reference_frequency;
% params_tx.power = 1.0e-3;
% params_tx.linewidth = 0;
% params_tx.bit_pattern = bit_pattern;
% params_tx.rise_time = 1/symbol_rate/4;
% sig = tx_ook_nrz(params_tx);
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
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

nrz_data_sig = elec_pulse_sequence_nrz(params.bit_pattern,params.rise_time);
% Generate NRZ data stream

sig = opt_laser_cw(params);
% CW laser

vpi = 1.0;
driving_signal_1 = vpi/2*(nrz_data_sig - 0.5);
driving_signal_2 = -vpi/2*(nrz_data_sig - 0.5);
bias_1 = 1.5*vpi;
bias_2 = 0;
split_in = 0.5;
split_out = 0.5;
loss = 0;
sig = mod_mzm(sig,driving_signal_1,driving_signal_2,bias_1,bias_2,vpi,split_in,split_out,loss);
% Data modulator

end