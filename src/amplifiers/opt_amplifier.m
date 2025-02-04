function sig = opt_amplifier(sig,params)
% Basic system model of an optical amplifier
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function simulates the operation of an ideal optical amplifier with
% wavelength independent gain and noise figure and where eventual
% saturation is governed by the signal average power.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_optamp.mode = 'power';%'gain';%'saturation';
% params_optamp.pol = 'x';%'y';%'both';
% params_optamp.gain = 20;
% params_optamp.output_power = 10;
% params_optamp.noise_figure = 4;
% params_optamp.add_noise = 1;
% sig = opt_amplifier(sig,params_optamp);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
%
% params            amplifier parameters [structure]
%
%                       params.mode      
%                           operation mode of the amplifier [string]

%                           params.mode = 'gain'
%                               The amplifier is gain-controlled. 
%                               A maximum allowed output power value may 
%                               clamp the gain.

%                           params.mode = 'power'
%                               The amplifier is output power controlled.

%                           params.mode = 'saturation'
%                               The amplifier gain is specified according 
%                               to a small signal gain and a saturated 
%                               output power value. 
%                               The saturation occurs on an average power 
%                               basis and therefore the module is suited at
%                               emulating the behaviour of EDFA-type 
%                               amplifiers rather than SOAs in this mode.                     
%
%                       params.output_power
%                           amplifier output power, in dBm [real scalar]
%
%                           Corresponds to the total output power in the 
%                           'power' mode, to the output saturation power 
%                           in the 'saturation' mode, or to the maximum 
%                           allowed output power in the 'gain' mode. 
%                           In this latter case, it should be ensured the 
%                           value is high enough, otherwise the amplifier 
%                           will indeed have a saturation behaviour.
%   
%                       params.gain
%                           amplifier gain, in dB [real scalar]
%
%                           Corresponds to the amplifier gain in the 'gain'
%                           mode, or to the small signal gain in the 
%                           'saturation' mode. 
%                           Not used in the 'power' mode.
%
%                       params.noise_figure
%                           amplifier noise figure, in dB [real scalar]
%
%                           If params.noise_figure = 0 dB, no noise is 
%                           added by the amplifier.
%
%                       params.add_noise
%                           specifies whether noise is added to the signal 
%                           or not [0/1]
%
%                           params.add_noise = 0;
%                           params.add_noise = 1;
%
%                       params.pol
%                           specifies whether noise is added to one or two
%                           polarisations [string]
%
%                           params.pol = 'x';
%                           params.pol = 'y';
%                           params.pol = 'both';
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output optical signal [optical signal structure]
%                       This is the amplified + noisy output signal.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% CONSTANT              essential physical constants [structure]
%
% dt                    time samples separation, in s [real scalar]
%
% reference_frequency   reference frequency, in Hz [real scalar]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global CONSTANT
global reference_frequency 
global dt


nsamples = length(sig.x);
% Number of samples in the input optical signal

input_power = char_opt_average_power(sig);
% Calculate average input power, in W

% Calculation of the gain, depending on the amplifier control mode.
if strcmp(params.mode,'gain') 
    % Gain mode
    
    gain_linear = 10^(params.gain/10);
    % Calculate gain, in linear units  
    output_power_max = 10^(params.output_power/10)*1.0e-3;
    % Maximum output power, in W
    output_power = input_power*gain_linear;
    % Estimate the output power, in W    
    if output_power > output_power_max
        gain_linear = output_power_max/input_power;
    end
    % Check if the output power is larger than the maximum allowed output 
    % power and clamp the gain if this is the case    
    
elseif strcmp(params.mode,'power')
    % Output power mode
    
    output_power = 10^(params.output_power/10)*1.0e-3;
    % Calculate expected total signal output power, in W 
    gain_linear = output_power/input_power;
    % Calculate the corresponding gain, in linear units 
    
elseif strcmp(params.mode,'saturation')     
    % Saturation mode
    
    saturated_output_power = 10^(params.output_power/10)*1.0e-3;
    % Calculate saturated output power, in W
    g0 = 10^(params.gain/10);
    % Calculate linear small signal gain, in linear units
    ps = saturated_output_power*(g0 - 2)/(g0*log(2));
    % Calculate the intrinsic saturation power, in W
    
    func = @(g) [g - g0*exp(-(g - 1)*input_power/ps)];
    % Define implicit function that needs to be solved to find the large
    % signal gain
    
    gain_linear = fzero(func,g0);
    % Solve for the large signal gain, in linear units 
    
end
% Now the gain has been calculated

sig.x = sig.x * sqrt(gain_linear);
sig.y = sig.y * sqrt(gain_linear);
% Multiply the input field by the gain. The amplifier is assumed to have no
% polarization-dependent gain.

output_power_signal_calculated = char_opt_average_power(sig);
% Calculate average signal output power, in W


if params.add_noise == 1
    % Once the gain has been calculated, it is the turn of the noise
    nsp = 0.5*(10^(params.noise_figure/10)*gain_linear - 1)/(gain_linear - 1);
    % Spontaneous emission factor
    % Eq. (2.118) on page 100 in E. Desurvire, "Erbium doped fiber amplifiers, 
    % Principles and applications, Wiley, 2002.
    noise_spectral_density = nsp*(gain_linear - 1)*CONSTANT.h*reference_frequency;
    % Noise spectral density per polarization
    noise_bandwidth = 1/dt;
    % Calculate the noise bandwidth. It is taken equal to the sample rate, 
    % or the inverse of the time interval between consecutive samples.
    ase_variance = noise_spectral_density*noise_bandwidth/2;
    % ASE noise variance, per quadrature and per polarisation
    ase_x = sqrt(ase_variance)*randn(1,nsamples) + 1i*sqrt(ase_variance)*randn(1,nsamples);
    ase_y = sqrt(ase_variance)*randn(1,nsamples) + 1i*sqrt(ase_variance)*randn(1,nsamples);
    % Create complex Gaussian noise for both polarisations.

    if strcmp(params.pol,'x')
        % Noise is added only to the -x polarisation. This would correspond
        % to having a polariser parallel to -x placed after the optical
        % amplifier.
        sig.x = sig.x + ase_x;    
        
    elseif strcmp(params.pol,'y')
        % Noise is added only to the -y polarisation. This would correspond
        % to having a polariser parallel to -y placed after the optical
        % amplifier.
        sig.y = sig.y + ase_y; 
        
    elseif strcmp(params.pol,'both');
        % Noise is added to both polarisations.
        sig.x = sig.x + ase_x;
        sig.y = sig.y + ase_y;   
        
    end
    
    
end
% end of params.add_noise = 1 mode

end