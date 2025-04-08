function eye = calc_eye(sig,nsamples_per_symbol)
% Quick eye diagram calculation
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the eye diagram of a real signal. The eye is
% represented over 1 symbol duration only.
% It is assumed that the input signal vector contains an integer number of
% symbols, which should normally be guaranteed by the way signals are
% generated.
% The plotting functionality is completely left out of the function.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% eye = calc_eye(sig,nsamples_per_symbol);
% figure('Name','eye diagram')
% plot(time_array(1:nsamples_per_symbol),eye,'b');
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig                   input signal [real vector]
%
% nsamples_per_symbol   number of samples per symbol in the input signal
%                           [integer scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% eye                   eye diagram [real matrix]
%
%                           The eye diagram is a matrix with:
%                           - nsamples_per_symbol lines
%                           - nsymbols columns
%                           i.e. each column in the matrix represents a
%                           trace in the eye diagram.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

sig = sig(:).';
% Ensure the signal is a line vector

nsamples = length(sig);
% Number of samples in the signal to display

nsymbols = nsamples/nsamples_per_symbol;
% We take for granted this is an integer. Should be the case from the way
% we generate the signal.

eye = reshape(sig,nsamples_per_symbol,nsymbols);
% Eye diagram

end