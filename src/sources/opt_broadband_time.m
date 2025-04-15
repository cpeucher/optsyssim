function sig = opt_broadband_time(power)
% Optical broadband source (time domain generation)
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates a white broadband signal with specified power 
% per frequency component. It simply generates an impulse in the time 
% domain.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% power_dbm = 0;
% sig = opt_broadband_time(power_dbm); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% power             power of each frequency component, in dBm
%                       [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output optical signal [optical signal structure]
%
%                       The power as a function of frequency is uniform. 
%                       The signal consists of a single impulse in the 
%                       time domain. 
%                       The signal is polarised along the '-x' direction.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% time_array        time samples, in s [real vector]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global time_array


nsamples = length(time_array);
% Number of frequency samples

sig = struct;
% Define the output signal as a structure

sig.x = zeros(1,nsamples);
sig.x(1) = sqrt(1.0e-3*10^(power/10))*nsamples;
sig.y = zeros(1,nsamples);
% Create the output optical signal

end