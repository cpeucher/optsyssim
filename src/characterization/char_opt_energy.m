function ep = char_opt_energy(sig,time_interval)
% Energy of an optical signal in a specified time interval
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the energy of an optical signal within a
% specified time interval.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% ep = char_opt_energy(sig,[time_array(1),time_array(end)]);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
%
% time_interval     limits of the interval over which the pulse power 
%                       will be integrated, in s [real vector]
%
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% ep                energy of the pulse over the specified interval, in J
%                       [real scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% time_array            time samples, in s [real vector]
%
% dt                    time samples separation, in s [real scalar]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


global time_array
global dt


time_interval = sort(time_interval);
% Ensure the elements in time_interval are in increasing order

power = abs(sig.x).^2 + abs(sig.y).^2;
% Calculate signal power

istart = find(time_array >= time_interval(1), 1);
% Find index corresponding to the start of integration range
istop = find(time_array <= time_interval(2), 1, 'last' );
% Find index corresponding to the end of the integration range

ep = num_int1d_simpson(power(istart:istop),dt);
% Integrate the signal power over the integration range.

end