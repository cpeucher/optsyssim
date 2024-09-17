function sig = tx_dpsk_nrz_pm(params)
% NRZ-DPSK transmitter using a phase modulator
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates an ideal NRZ-DPSK signal using a phase modulator.
% The transmitter consists of a CW laser followed by a phase modulator.
% Note that differential encoding of the data to be modulated is not
% performed within this module and should therefore be performed, whenever
% needed, before calling the function.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_tx.emission_frequency = reference_frequency;
% params_tx.power = 1.0e-3;
% params_tx.linewidth = 0;
% params_tx.bit_pattern = logical_differential_encoder_binary(bit_pattern);
% params_tx.rise_time = 1/symbol_rate/4;
% sig = tx_dpsk_nrz_pm(params_tx); 
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
% CW laser signal.

sig = mod_pm(sig,nrz_data_sig,1,0);
% Phase modulator.

end