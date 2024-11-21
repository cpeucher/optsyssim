function peak_power = char_opt_peak_power(sig,time_interval)
% Peak power of optical signal
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the peak power of an optical signal within a
% specified time interval.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% peak_power = char_opt_peak_power(sig,time_interval);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
%
% time_interval     2 element vector specifying the limits of the interval
%                       over which the pulse peak power will be searched,
%                       in s [real vector]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% peak_power        peak power of the signal in the specified interval, 
%                       in W [real scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% time_array        time samples, in s [real vector]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global time_array

power = abs(sig.x).^2+abs(sig.y).^2;
% Calculate pulse power

time_interval = sort(time_interval);
% Ensure the elements in time_interval are in increasing order

i_start = min(find(time_array >= time_interval(1)));
% Find index corresponding to the start of integration range
i_stop = max(find(time_array <= time_interval(2)));
% Find index corresponding to the end of the integration range

peak_power = max(power(i_start:i_stop));
% Find peak power in the range of interest

end