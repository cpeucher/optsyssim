function sig = opt_source_pulse(symbs,params)
% Optical pulse source
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates an optical pulse train according to a specified
% binary sequence. The generated signal is polarised along -x.
% The pulse train generation method consists in convolving a Dirac
% comb with proper complex amplitudes with the desired pulse shape.
% Therefore complex modulation of the optical pulse train is possible.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_pulse_train.type = 'gaussian';%'sech';
% params_pulse_train.order = 1;
% params_pulse_train.emission_frequency = reference_frequency + 100e9;
% params_pulse_train.peak_power = 1e-3;
% params_pulse_train.fwhm = 50e-12;
% params_pulse_train.chirp = 5;
% sig = opt_source_pulse(seq,params_pulse_train); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% symbs             complex symbols sequence according to which the pulse
%                       train will be modulated [complex vector]
%
% params            pulse train parameters [structure]
%
%                       params.type
%                           type of pulse [string]
%                           params.type = 'gaussian';
%                           params.type = 'sech';
%
%                       params.order
%                           Gaussian order [integer scalar]
%
%                       params.peak_power
%                           pulse peak power, in W [real scalar]
%
%                       params.fwhm
%                           pulse full-width at half maximum duration, in s 
%                           [real scalar]
%
%                       params.chirp
%                           pulse linear chirp parameter, no unit 
%                           [real scalar]
%
%
%                       params.emission_frequency
%                           desired signal centre frequency, in Hz 
%                           [real scalar]
%                           The frequency may be adjusted so that it falls
%                           on one sample of the frequency grid of the
%                           simulation.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               modulated optical signal [optical signal structure]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% time_array            time samples, in s [real vector]
%
% reference_frequency   reference frequency, in Hz [real scalar]
%
% df                    frequency samples separation, in Hz [real scalar]
%
% -------------------------------------------------------------------------
% REMARKS:
% -------------------------------------------------------------------------
% 
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

global time_array 
global reference_frequency
global df 


nsamples = length(time_array);
% Number of time samples.
nsymbols = length(symbs);
% Number of symbols.

if rem(nsamples,nsymbols) == 0
    % Check that the number of samples is an integer multiple of the number
    % of symbols.
    samples_per_symbol = nsamples/nsymbols;
    % Number of samples per symbol.
else
    error('elec_modulator: the number of samples should be an integer multiple of the number of symbols.');
end


position = time_array(length(time_array)/2);
% Position of the isolated pulses that will be used in the convolution with
% the Dirac comb.

% First we generate the isolated optical pulses shapes.
switch params.type    
        
    case 'sech'
        % sech pulse.
        
        pulse = opt_pulse_sech(time_array,params.peak_power,position,params.fwhm,params.chirp);
        
    case 'gaussian'
        % Gaussian pulse.
              
        pulse = opt_pulse_gaussian(time_array,params.order,params.peak_power,position,params.fwhm,params.chirp);
        
    otherwise
        
        error('opt_source_pulse: pulse shape not defined.');
        
end

seq = zeros(1,nsamples);
% Create an empty array for the upsampled sequence.
seq(samples_per_symbol/2:samples_per_symbol:end) = symbs;
% Upsample the input data symbols with 1 sample per symbol and zeros
% in-between.

pulse = ifft(fft(seq).*fft(pulse));
% Modulated signal. Convolution of the Dirac comb and the pulse shape.

pulse = circshift(pulse,[0 -nsamples/2+1]);
% Compensate delay.

% Until know, the generated signal is represented by a complex vector.
% We need to convert it to the optical signal format, and translate the
% envelope to the desired frequency.


fprintf(1,'\n%s\n%s%f%s\n','opt_source_pulse:','Desired emission frequency: ',params.emission_frequency/1.0e12,' THz');
% Display the desired emission frequency.

actual_emission_frequency = reference_frequency + round((params.emission_frequency - reference_frequency)/df)*df;
% First we ensure that the emission frequency falls onto the frequency
% grid.
% Not really necessary for a pulsed signal with broad spectrum. But anyway,
% done for the sake of good order...



fprintf(1,'%s%f%s\n','Actual emission frequency: ',actual_emission_frequency/1.0e12,' THz');
fprintf(1,'%s%f%s\n','Delta f: ',(actual_emission_frequency - params.emission_frequency)/1.0e9,' GHz');
% Ee notify the the user of this minor change.


sig = [];
sig.x = pulse.*exp(1i*2*pi*(actual_emission_frequency - reference_frequency)*time_array);
sig.y = zeros(1,nsamples);
% Generate optical signal.
    
end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------