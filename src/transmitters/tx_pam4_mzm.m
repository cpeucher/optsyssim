function sig = tx_pam4_mzm(params)
% PAM4 transmitter using a Mach-Zehnder modulator (chirp-free, infinite ER)
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates an optical PAM4 signal from two binary sequences
% using a dual-drive Mach-Zehnder modulator.
% The generation of the PAM4 signal is limited to chirp-free (push-pull) 
% modulation with infinite extinction ratio (utilising the full-swing of 
% the Mach-Zehnder modulator).
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_tx.emission_frequency = reference_frequency;
% params_tx.linewidth = 0;
% params_tx.power = 1.0e-3;
% params_tx.bit_pattern_1 = bit_pattern_1;
% params_tx.bit_pattern_2 = bit_pattern_2;
% params_tx.rise_time = 1/symbol_rate/4;
% sig = tx_pam4_mzm(params); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% params            transmitter parameters [structure]
%
%                       params.emission_frequency
%                           desired centre frequency of the channel, in Hz
%                           [real scalar]
%                           It will be adapted so that it coincides 
%                           with a sample of the discrete frequency grid.
%
%                       params.linewidth
%                           laser linewidth, in Hz [real scalar]
%
%                       params.power
%                           laser power, in W [real scalar]
%
%                       params.bit_pattern_1
%                           binary data pattern [binary vector]
%
%                       params.bit_pattern_2
%                           binary data pattern [binary vector]
%
%                       params.symbol_rate
%                           symbol rate of the signal to be generated, 
%                           in baud [real scalar]
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

nrz_data_sig_1 = elec_pulse_sequence_nrz(params.bit_pattern_1,params.rise_time);
nrz_data_sig_2 = elec_pulse_sequence_nrz(params.bit_pattern_2,params.rise_time);
% Generate NRZ data streams

sig = opt_laser_cw(params);
% CW laser signal

alpha = 0.5 + asin(1/3)/pi;
rho = (0.5 - asin(1/3)/pi)/(0.5 + asin(1/3)/pi);
% Parameters for ensuring full-swing modulation with equally-spaced power
% levels. See below.

pam4_data_sig = alpha*(nrz_data_sig_1 + rho*nrz_data_sig_2);
% Modulating signal within [0 1]

vpi = 1.0;
vpp = vpi/2;
driving_signal_1 = vpp*(pam4_data_sig-0.5);
driving_signal_2 = -vpp*(pam4_data_sig-0.5);
bias_1 = 1.5*vpi;
bias_2 = 0;
sig = mod_mzm(sig,driving_signal_1,driving_signal_2,bias_1,bias_2,vpi,0.5,0.5,0);
% Mach-Zehnder modulator

end