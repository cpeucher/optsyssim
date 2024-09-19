function sig = elec_rise_time(sig,rise_time)
% Apply rise time to a rectangular electrical signal
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function applies a rise time to an electrical rectangular signal.
% The rise time is obtained by filtering the original signal with a 1st
% order Gaussian electrical filter.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% sig = elec_rise_time(sig,rise_time);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input electrical signal [real vector]
%
% rise_time         rise time to impose to the input signal, in s 
%                       [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output electrical signal with rise time [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

rise_time_filter.type = 'gaussian';
rise_time_filter.order = 1;
rise_time_filter.f3dB = 0.3321/rise_time;
% Relation between cut-off frequency and rise time for a 1st order Gaussian
% electrical low-pass filter.       
    
sig = elec_elpf(sig,rise_time_filter);
% Apply rise time filter.


end