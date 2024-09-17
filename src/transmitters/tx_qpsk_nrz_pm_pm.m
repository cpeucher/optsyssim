function sig = tx_qpsk_nrz_pm_pm(params)
% NRZ-QPSK transmitter using two phase modulators
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates an ideal NRZ-QPSK signal using two concatenated
% phase modulators.
% The first phase modulaor produces a (0,pi) phase shift while the second 
% modulator produces a (0,pi/2) phase shift.
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
% sig = tx_qpsk_nrz_pm_pm(params_tx); 
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

nrz_data_sig_pm_1 = elec_pulse_sequence_nrz(params.bit_pattern_1,params.rise_time);
% NRZ data stream to be applied to the first PM
nrz_data_sig_pm_2 = elec_pulse_sequence_nrz(params.bit_pattern_2,params.rise_time);
% NRZ data stream to be applied to the second PM

sig = opt_laser_cw(params);
% CW laser signal

sig = mod_pm(sig,nrz_data_sig_pm_1,1,0);
% 1st phase phase modulator

sig = mod_pm(sig,nrz_data_sig_pm_2/2,1,0);
% 2nd phase phase modulator

end