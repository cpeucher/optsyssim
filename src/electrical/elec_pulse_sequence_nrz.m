function sig = elec_pulse_sequence_nrz(seq,rise_time)
% Electrical NRZ pulse sequence generation
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates a sequence of NRZ electrical pulses. The sequence 
% follows the binary data in the array seq. 
% The bit rate is determined from the length of the seq vector and the
% duration of the simulation time window. 
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% sig = elec_pulse_sequence_nrz(bit_pattern,rise_time);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% seq               input binary sequence [binary vector]
%
% rise_time         rise time of the pulses (10-90%) [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output electrical waveform [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% time_array        time samples, in s [real vector]
%
% dt                time samples separation, in s [real scalar]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global time_array
global dt


nsymbols = length(seq);
% Number of bits

nsamples = length(time_array);
% Number of samples

nsps = nsamples/nsymbols;
% Number of samples per bit

if mod(nsps,1) ~= 0 
    % Check if nsps is an integer
    error('elec_pulse_sequence_nrz: non-integer number of samples per symbol! Adjust the sequence length.');
end

fprintf('\n\n%s\n%s%3.3e%s\n','elec_pulse_sequence_nrz:','Bit rate of sequence is ',1/nsps/dt,' bit/s.')

sig = seq(ones(1,nsps),:);
sig = sig(:).';    
% Create ideal rectangular NRZ pulse sequence    

sig = elec_rise_time(sig,rise_time);
% Apply rise time filter
            
end 