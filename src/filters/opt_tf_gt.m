function tf = opt_tf_gt(freq,params)
% Transfer function of Gires-Tournois interferometer
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the transfer function of an ideal Gires-Tournois
% interferometer.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_gt.reflectivity = 0.9;
% params_gt.fsr = 100e9;
% params_gt.centre_frequency = 0;
% tf = opt_tf_gt(frequency_array,params_gt);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% params            parameters of the GT filter [structure]
%
%                       params.reflectivity
%                           intensity reflection coefficient of the input
%                           interface [real scalar]
%
%                       params.fsr
%                           free spectral range, in Hz [real scalar]
%
%                       params.centre_frequency
%                           centre frequency of one of the resonances,
%                           in Hz [real scalar]
%
%                           This value can be specified as absolute or
%                           relative (with respect to the centre frequency
%                           of the simulation bandwidth, i.e. the global
%                           parameter reference_frequency), depending on
%                           whether freq is specified as absolute or
%                           relative frequency.
%
% freq              frequencies at which the transfer function 
%                       will be calculated, in Hz [real vector]
%
%                       This can be specified as either absolute or
%                       relative frequencies, but in a way that is
%                       consistent with the params.centre_frequency
%                       parameter.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% tf                transfer function of the GT interferometer
%                       [complex vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

phi = 2*pi*(freq - params.centre_frequency)/params.fsr;
% Phase shift term.

tf = (sqrt(params.reflectivity) - exp(-1i*phi))./(1 - sqrt(params.reflectivity)*exp(-1i*phi));
% Transfer function.  

end