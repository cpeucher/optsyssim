function [sig, dclevel] = elec_dc_block(sig)
% Electrical DC block
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This functions removes the dc content of an electrical signal.
% It is based on the calculation of the mean value of the signal. It is
% therefore not robust when the signal is imbalanced (for instance when the
% numbers of marks and spaces differ significantly).
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% [sig, dclevel] = elec_dc_block(sig); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input electrical signal [real vector]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output electrical signal [real vector]
%
% dclevel           estimated dc level [real scalar]
% 
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

nsamples = length(sig);
% Number of samples in the signal.

fseries = fft(sig)/nsamples;
% Calculates the spectrum of the signal

dclevel = fseries(1);

sig = sig - dclevel;
% Removes the dc content from the signal

end