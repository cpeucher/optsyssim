function sig = elec_coder_manchester(seq,nsamples_per_symbol)
% Manchester encoder
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements a Manchester electrical encoder.
% The generated signal takes values between 0 and 1 (unbalanced).
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% sig = elec_coder_manchester(bit_pattern,nsamples_per_symbol); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% seq                   data sequence to modulate [binary vector]
%
% nsamples_per_symbol   number of samples per bit of the modulated signal
%                           [integer scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig                   electrical modulated signal [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

nrz_data_sig = elec_coder_nrz(seq,nsamples_per_symbol);
% Generate NRZ sequence.

sig_clock = [ones(1,nsamples_per_symbol/2) zeros(1,nsamples_per_symbol/2)];
% One clock period.
sig_clock = repmat(sig_clock,1,length(seq));
% Clock signal with 50% duty cycle.

sig = xor(nrz_data_sig,sig_clock);
% The Manchester signal is the result of the XOR operation between the
% clock and the NRZ signal.

end