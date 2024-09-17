function sig = tx_dpsk_rz67_mzm(params)
% 67% RZ-DPSK transmitter using a Mach-Zehnder modulator
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates an ideal RZ67-DPSK signal using a Mach-Zehnder
% modulator.
% The transmitter consists of a CW laser followed by a Mach-Zehnder data 
% modulator. The generated signal is ideal in the sense that the extinction
% ratios of the data modulator is infinite, and chirp-free modulation is 
% applied.
% The data modulator is biased at minimum point and driven with a 
% peak-to-peak voltage difference between the 2 driving signals equal to 
% 2Vpi.
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
% sig = tx_dpsk_rx67_mzm(params_tx); 
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
% dt                time samples separation, in s [real scalar]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global dt


nrz_data_sig = elec_pulse_sequence_nrz(params.bit_pattern,params.rise_time);
% Generate NRZ data stream

bit_rate = length(params.bit_pattern)/length(nrz_data_sig)/dt;
% Signal bit rate

params_rf.frequency = bit_rate/2;
params_rf.phase = pi/2;
params_rf.vpp = 1.0;
params_rf.vdc = 0;
clock_sig = elec_sinusoidal(params_rf);
% Generate electrical clock signal

sig = opt_laser_cw(params);
% CW laser signal

vpi = 1.0;
driving_signal_1 = vpi*(nrz_data_sig - 0.5);
driving_signal_2 = -vpi*(nrz_data_sig - 0.5);
bias_1 = vpi;
bias_2 = 0;
sig = mod_mzm(sig,driving_signal_1,driving_signal_2,bias_1,bias_2,vpi,0.5,0.5,0);
% Data modulator

bias_1 = vpi;
bias_2 = 0;
driving_signal_1 = vpi*clock_sig;
driving_signal_2 = -vpi*clock_sig;
sig = mod_mzm(sig,driving_signal_1,driving_signal_2,bias_1,bias_2,vpi,0.5,0.5,0);
% Pulse carver

end