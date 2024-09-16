function sig = opt_combiner_y_junction(sig1,sig2)
% Optical Y-junction combiner with 50% combining power ratio
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function models an ideal 1x2 Y-junction combiner. The Y junction is
% assumed to be ideal with power combining ratio equal to 1/2. The state of
% polarisation is assumed to be preserved by the Y-junction.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% sig = opt_combiner_y_junction(sig1,sig2);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig1              input optical signal at port 1
%                       [optical signal structure]
%
% sig1              input optical signal at port 2
%                       [optical signal structure]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output optical signal [optical signal structure]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

sig.x = sqrt(1/2)*sig1.x + sqrt(1/2)*sig2.x;
sig.y = sqrt(1/2)*sig1.y + sqrt(1/2)*sig2.y;

end