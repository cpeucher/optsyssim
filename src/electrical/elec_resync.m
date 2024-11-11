function [sig, delay_samples] = elec_resync(sig,sig_ref)
% Resynchronise electrical signal to a reference signal
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function resynchronises an electrical signal with respect to a 
% reference electrical signal.
% This can be useful in order to evaluate performance metrics such as BER
% etc, represent eye diagram etc.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% [sig, delay_samples] = elec_resync(sig,sig_ref);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               electrical signal to retime [real vector]
%
% sig_ref           reference electrical signal [real vector]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               retimed electrical [real vector]
%
% delay_samples     delay, expressed in terms of number of samples, between
%                       the original and the retimed signal (or the
%                       reference signal) [integer scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

[~,delay_samples] = max(xcorr(sig,sig_ref));
% Calculate the delay shift that maximises the cross-correlation between 
% the signal and the reference signal

sig = circshift(sig,[0 -delay_samples]);
% Retime the original signal

end