function tf = opt_tf_fp(freq,params)
% Transfer function of Fabry-Perot filter
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the transfer function of an ideal Fabry-Perot 
% filter.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_fp.mode = 'transmission';%'reflection';
% params_fp.finesse = 100;
% params_fp.fsr = 100e9;
% params_fp.centre_frequency = 0;
% tf = opt_tf_fp(frequency_array,params_fp);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% params            parameters of the Fabry-Perot filter [structure]
%
%                       params.mode
%                           specifies whether the Fabry-Perot filter is
%                           used in transmission or in reflection [string]
%
%                           params.mode = 'transmission';
%                           params.mode = 'reflection';
%
%                       params.finesse
%                           finesse of the Fabry-Perot filter [real scalar]
%
%                       params.fsr
%                           free spectral of the Fabry-Perot filter, in Hz
%                           [real scalar]
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
% tf                transfer function of the Fabry-Perot filter
%                       [complex vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

a = params.finesse*params.finesse;
b = -(2*a + pi*pi);
c = a;
R = (-b+sqrt(b*b-4*a*c))/(2*a);
% Calculate the energy reflectivity from the desired finesse.

T = 1 - R;
% Energy transmittivity.

phi = 2*pi*(freq - params.centre_frequency)/params.fsr;
% Phase shift term.


if strcmp(params.mode,'transmission') 
    
    tf = T./(1 - R*exp(-1i*phi));    
    
elseif strcmp(params.mode,'reflection')
    
    tf = sqrt(R)*(1 - exp(-1i*phi))./(1 - R*exp(-1i*phi));
    
else
    disp('opt_tf_fp: mode should be either transmission or reflection.')    
end

end