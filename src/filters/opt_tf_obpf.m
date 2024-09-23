function tf = opt_tf_obpf(params,freq)
% Transfer functions of some standard optical bandpass filters
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the transfer functions of optical bandpass
% filters of various standard types.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_obpf.type = 'gaussian';%'none','rectangular_ideal','rectangular';
% params_obpf.centre_frequency = 0;
% params_obpf.bandwidth = 40e9;
% params_obpf.order = 4;
% params_obpf.attenuation_in_band = -1;
% params_obpf.attenuation_out_band = -20;
% tf = opt_tf_obpf(params_obpf,frequency_array);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% params            optical filter parameters [structure]
%
%                       params.type
%                           type of optical filter [string]
%
%                           params.type = 'none'
%                               no filter (all-pass)
%                               Allows to have the filter function in a 
%                               script and turn the filter 'off' if
%                               necessary.
%
%                           params.type = 'gaussian'
%                               (super-)Gaussian filter transfer function
%                               This is a Gaussian function as a function
%                               of frequency (not wavelength)
%
%                           params.type = 'rectangular_ideal'
%                               rectangular filter with transmission 1 in
%                               the passband and 0 out-of-band
%
%                           params.type = 'rectangular'
%                               rectangular filter where the in-band and 
%                               out-of-band transmission can be specified 
%                               through the 
%                                   params.attenuation_in_band and 
%                                   params.attenuation_out_band
%                               parameters, respectively
%                               Can be used for bandpass or bandstop 
%                               (notch) filters.
%
%                       params.order
%                           Gaussian order [integer scalar]
%
%                           For params.type = 'gaussian' only
%
%                       params.bandwidth 
%                           FWHM bandwidth (in frequency, defined at half
%                           the maximum value of the power transfer 
%                           function), in Hz [real scalar]
%
%                       params.centre_frequency
%                           centre frequency of the passband, in Hz
%                           [real scalar]
%                           Can be expressed in terms of absolute frequency
%                           (e.g. 193.1e12 Hz), or in terms of relative
%                           frequency (e.g. -40e9 Hz), depending on whether
%                           the input vector freq is expressed in absolute
%                           or relative frequency.
%                           In case freq is the default sampling frequency
%                           array, i.e. freq = 'frequency_array', then 
%                           params.centre_frequency = 0 means that the 
%                           filter is tuned to the centre frequency of the 
%                           array, i.e. to reference_frequency.    
%
%                       params.attenuation_in_band
%                           in-band attenuation, in dB (negative number for 
%                           loss) [real scalar]
%                           For params.type = 'rectangular' only.    
%
%                       params.attenuation_out_band
%                           out-of-band attenuation, in dB (negative number 
%                           for loss) [real scalar]
%                           For params.type = 'rectangular' only.  
%
% freq              frequency values at which the transfer function is 
%                       calculated, in Hz
%
%                       Can be relative of absolute frequencies depending
%                       on whether params.centre_frequency is a relative or
%                       absolute frequency, when this parameter is
%                       relevant, or whether some input frequency values
%                       are provided (for instance for params.type =
%                       'custom').
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% tf                values of the transfer function that applies to the 
%                       field [complex vector]
%                       Square for power transfer function.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%--------------------------------------------------------------------------
% ALL-PASS
%--------------------------------------------------------------------------
if strcmp(params.type,'none')
    % No filtering is applied.
    
    tf = ones(1,length(freq));    
    
%--------------------------------------------------------------------------
% GAUSSIAN
%--------------------------------------------------------------------------
elseif strcmp(params.type,'gaussian')
    % Gaussian optical bandpass filter.
       
    a = params.bandwidth/2/(log(2))^(1/(2*params.order));
    % Conversion of FWHM bandwidth to half-bandwidth at 1/e.
    tf = exp(-0.5*((freq - params.centre_frequency)/a).^(2*params.order));    
    
%--------------------------------------------------------------------------
% RECTANGULAR
%--------------------------------------------------------------------------    
elseif strcmp(params.type,'rectangular_ideal')   
    % Rectangular optical bandpass filter.
    
    tf = abs(freq - params.centre_frequency) <= params.bandwidth/2;

%--------------------------------------------------------------------------
% RECTANGULAR-NOTCH
%--------------------------------------------------------------------------
elseif strcmp(params.type,'rectangular')
   
    in_band = abs(freq - params.centre_frequency) <= params.bandwidth/2;
    out_band = abs(freq - params.centre_frequency) > params.bandwidth/2;
    
    tf = in_band*10^(params.attenuation_in_band/20) + out_band*10^(params.attenuation_out_band/20);
    
%--------------------------------------------------------------------------
% OTHERWISE...
%--------------------------------------------------------------------------
else
    error('opt_tf_obpf: optical filter type not implemented.');
end


end