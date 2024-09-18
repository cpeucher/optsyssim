function sig = tx_qpsk_nrz_iq(params)
% NRZ-QPSK transmitter using optical IQ modulator
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates an ideal NRZ-QPSK signal using an I/Q modulator.
% The transmitter consists of a CW laser followed by an I/Q data modulator. 
% The generated signal is ideal in the sense that the extinction
% ratios of the data modulator is infinite and chirp-free modulation is 
% applied.
% The data modulators are biased at minimum point and driven with a 
% peak-to-peak voltage difference between the 2 driving signals equal to 
% 2Vpi.
% Note that, for DQPSK, differential encoding of the data to be modulated
% is not performed within this module and should therefore be performed, 
% whenever needed, before calling the function.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% u = round(rand(1,nsymbols));
% v = round(rand(1,nsymbols));
% [params_tx.bit_pattern_1,params_tx.bit_pattern_2,~,~,~] = logical_differential_encoder_dqpsk('parallel',u,v);
% params_tx.emission_frequency = reference_frequency;
% params_tx.linewidth = 0;
% params_tx.power = 1.0e-3;
% params_tx.rise_time = 1/symbol_rate/4;
% sig = tx_qpsk_nrz_iq(params_tx); 
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

nrz_data_sig_i = elec_pulse_sequence_nrz(params.bit_pattern_1,params.rise_time);
nrz_data_sig_q = elec_pulse_sequence_nrz(params.bit_pattern_2,params.rise_time);
% Generate NRZ data streams for I and Q quadratures

sig = opt_laser_cw(params);
% CW laser signal

vpi = 1.0;
driving_signal_i1 = vpi*(nrz_data_sig_i - 0.5);
driving_signal_i2 = -vpi*(nrz_data_sig_i - 0.5);
driving_signal_q1 = vpi*(nrz_data_sig_q - 0.5);
driving_signal_q2 = -vpi*(nrz_data_sig_q - 0.5);
bias_1 = vpi;
bias_2 = 0;


[sig_1,sig_2] = opt_splitter_y_junction(sig);
% Input power splitter
sig_1 = mod_mzm(sig_1,driving_signal_i1,driving_signal_i2,bias_1,bias_2,vpi,0.5,0.5,0);
% Data modulator for I quadrature
sig_2 = mod_mzm(sig_2,driving_signal_q1,driving_signal_q2,bias_1,bias_2,vpi,0.5,0.5,0);
% Data modulator for Q quadrature
sig_2 = opt_phase_shift(sig_2,pi/2);
% pi/2 phase shift applied to the Q quadrature
sig = opt_combiner_y_junction(sig_1,sig_2);
% Output power combiner

end