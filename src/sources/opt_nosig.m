function sig = opt_nosig()
% Dummy optical signal generation
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates an optical signal with zero power. This module is
% to be used as an empty input to modules requiring multiple inputs, such 
% as e.g. couplers.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% sig = opt_nosig();
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% None             
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output optical signal with null power
%                       [optical signal structure]
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
% Number of samples in the optical null signal.
sig = struct;
% Create signal structure.
sig.x = zeros(1,nsamples);
sig.y = zeros(1,nsamples);
% The field along the two polarisation states is null.

end