function sig = elec_coder_nrz(seq,nsamples_per_symbol)
% Electrical NRZ encoder without rise time
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function creates a rectangular NRZ signal from a binary sequence. It
% is essentially contained in the elec_pulse_sequence.m routine but is
% reproduced here for simplicity.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% sig = elec_coder_nrz(bit_pattern,nsamples_per_symbol); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% seq                   binary sequence to be modulated [binary vector]
%
% nsamples_per_symbol   number of samples per bit in the generated signal
%                           [integer scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig                   rectangular (no rise time) electrical NRZ signal 
%                           encoded by the sequence seq [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

sig = seq(ones(1,nsamples_per_symbol),:);

% One could also have done:
% sig = kron(seq,ones(1,samples_per_bit));
% which seems to be slower than the finally adopted solution.

sig = sig(:).';
% Ensure the signal is a saved as a line vector.

end