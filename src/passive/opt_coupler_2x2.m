function [sig21,sig22] = opt_coupler_2x2(sig11,sig12,mode,coupling)
% 2x2 optical directional coupler
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function models an ideal optical directional 2x2 coupler. The pi/2
% relative phase shift between the 2 output ports is taken into account.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% [sig21,sig22] = opt_coupler_2x2(sig11,sig12,'lin',0.5); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig11             input optical signal, upper port
%                       [optical signal structure]
%
% sig12             input optical signal, lower port
%                       [optical signal structure]
%
% mode              way the coupling ratio is specified [string]
%   
%                       mode = 'lin' means that the coupling factor is 
%                           specified as 0 <= coupling <= 1
%                           coupling = 0.5 for a 3 dB coupler   (50%-50%)
%                           coupling = 0.1 for a 10 dB coupler  (90%-10%)
%                           coupling = 0.01 for a 20 dB coupler (99%-1%)
%
%                       mode = 'log' means that the coupling factor is 
%                           specified in dB
%                           coupling = 3.0103 for a 3 dB coupler  (50%-50%)
%                           coupling = 10 for a 10 dB coupler     (90%-10%)
%                           coupling = 20 for a 20 dB coupler     (99%-1%)
%                           coupling = coupler coupling ratio
%
% coupling          coupling ratio, in linear or logarithmic units
%                       [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig21             output optical signal, upper port 
%                       [optical signal structure]
%
% sig22             output optical signal, lower port
%                       [optical signal structure]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

if strcmp(mode,'lin')
    alpha = 1 - coupling;
    
elseif strcmp(mode,'log')
    alpha = 1 - 10^(-coupling/10);
    
end    

sig21.x = sqrt(alpha)*sig11.x - 1i*sqrt(1 - alpha)*sig12.x;
sig21.y = sqrt(alpha)*sig11.y - 1i*sqrt(1 - alpha)*sig12.y;
sig22.x = -1i*sqrt(1 - alpha)*sig11.x + sqrt(alpha)*sig12.x;
sig22.y = -1i*sqrt(1 - alpha)*sig11.y + sqrt(alpha)*sig12.y;

end