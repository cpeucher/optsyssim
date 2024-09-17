function sig = tx_qpsk_nrz_mzm_pm(params)
% NRZ-QPSK transmitter using a Mach-Zehnder modulator followed by a phase modulator
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates an ideal NRZ-QPSK signal using a Mach-Zehnder
% modulator followed by a phase modulator.
% The Mach-Zehnder modulator is biased at minimum point and driven with a 
% peak-to-peak voltage difference between the 2 driving signals equal to 
% 2Vpi, resulting in a pi phase modulation. Its extinction ratio is 
% infinite. The modulation index of the phase modulator corresponds to 
% Vpi/2.
% Note that, for DQPSK, differential encoding of the data to be modulated
% is not performed within this module and should therefore be performed, 
% whenever needed, before calling the function.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% u = round(rand(1,nsymbols));
% v = round(rand(1,nsymbols));
% [params_tx.bit_pattern_1,params_tx.bit_pattern_2,~,~,~] = logical_differential_encoder_dqpsk('serial',u,v);
% params_tx.emission_frequency = reference_frequency;
% params_tx.linewidth = 0;
% params_tx.power = 1.0e-3;
% params_tx.rise_time = 1/symbol_rate/4;
% sig = tx_qpsk_nrz_mzm_pm(params_tx);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% params            transmitter parameters [structure]
%
%                   params.emission_frequency
%                       desired centre frequency of the channel, in Hz
%                       [real scalar]
%                       It will be adapted so that it coincides with a 
%                       sample of the discrete frequency grid.
%
%                   params.linewidth
%                       laser linewidth, in Hz [real scalar]
%
%                   params.power
%                       laser power, in W [real scalar]
%
%                   params.bit_pattern_1
%                       binary data pattern for the in-phase (I) component
%                       [binary vector]
%
%                   params.bit_pattern_2
%                       binary data pattern for the quadrature (Q) 
%                       component [binary vector]
%
%                   params.rise_time
%                       rise time of the electrical NRZ signal driving the 
%                       data modulator, in s [real scalar]
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

nrz_data_sig_mzm = elec_pulse_sequence_nrz(params.bit_pattern_1,params.rise_time);
% NRZ data stream to be applied to the MZM
nrz_data_sig_pm = elec_pulse_sequence_nrz(params.bit_pattern_2,params.rise_time);
% NRZ data stream to be applied to the PM

sig = opt_laser_cw(params);
% CW laser signal

vpi = 1.0;
driving_signal_1 = vpi*(nrz_data_sig_mzm - 0.5);
driving_signal_2 = -vpi*(nrz_data_sig_mzm - 0.5);
bias_1 = vpi;
Bias_2 = 0;
sig = mod_mzm(sig,driving_signal_1,driving_signal_2,bias_1,Bias_2,vpi,0.5,0.5,0);
% MZ modulator

sig = mod_pm(sig,nrz_data_sig_pm/2,vpi,0);
% Phase modulator

end