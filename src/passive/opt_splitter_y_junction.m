function [sig1,sig2] = opt_splitter_y_junction(sig)
% Y-junction splitter with 50% power splitting ratio
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function models an ideal 1x2 Y-junction splitter. The Y junction is
% assumed to be ideal with power splitting ratio equal to 1/2. The state of
% polarisation is assumed to be preserved by the Y-junction.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig1              optical signal at output port 1
%                       [optical signal structure]
%
% sig2              optical signal at output port 2
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

sig1.x = sqrt(1/2)*sig.x;
sig1.y = sqrt(1/2)*sig.y;
sig2.x = sqrt(1/2)*sig.x;
sig2.y = sqrt(1/2)*sig.y;

end