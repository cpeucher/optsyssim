function sig = opt_att(sig,params)
% Optical attenuator
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function attenuates an optical signal. Both polarisations are
% attenuated the same amount.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_att.mode = 'attenuation';%'power';
% params_att.target = 10;  % attenuation in dB or total output power in dBm
% sig = opt_att(sig,params_att);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
%
% params            attenuator parameters [structure]
%
%                       params.mode
%                           mode of the attenuator [string]
%
%                               params.mode = 'attenuation' means that the
%                                   attenuator attenuates the signal by the
%                                   amount specified by the parameter
%                                   target.
%
%                               params.mode = 'power' means that the 
%                                   attenuator will attenuate the signal so
%                                   that a total average power equal to 
%                                   target is obtained at the output.
%
%                       params.target
%                           target of the attenuator [real scalar]
%
%                           In 'attenuation' mode: attenuation, in dB. 
%                           For attenuation a positive number should be 
%                           provided. Otherwise the function will actually 
%                           amplify the signal.
%
%                           In 'power' mode: total output power, in dBm. 
%                           For attenuation a value smaller than the total 
%                           input power should be provided. Otherwise the
%                           function will actually amplify the signal.
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
% REMARKS:
% -------------------------------------------------------------------------
% The function may implement gain. No warning is issued. 
%
% -------------------------------------------------------------------------
% TO DO:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% CREDITS:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% AUTHOR:
% -------------------------------------------------------------------------
% Christophe Peucheret (christophe.peucheret@univ-rennes.fr)
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

if strcmp(params.mode,'attenuation')
    
    att = params.target;    
    
elseif strcmp(params.mode,'power')
    
    pin = char_opt_average_power(sig);
    % Input power to the attenuator, in W.
    pout = 10^(params.target/10)*1.0e-3;
    % Wished output power, in W
    att = -10*log10(pout/pin);
    % Corresponding attenuation, in dB. 
    % Should be a positive number for attenuation.
else
    error('opt_att: attenuation mode not implemented.');
end

sig.x = sig.x*10^(-att/20);
sig.y = sig.y*10^(-att/20);
% Attenuate the -x and -y polarisations of the signal field by the same
% amount.

end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------