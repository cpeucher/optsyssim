function sig = elec_modulator(symbs,params)
% Electrical modulator using convolution of impulses with pulse shape
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates electrical modulated waveforms using classic 
% pulse shapes such as raised-cosine, root raised-cosine, sinc, sech,
% Gaussian etc.
% The pulses can be modulated using complex modulation formats following
% proper mapping. The signal is generated by convolving the pulse shape
% with a Dirac comb of complex amplitude with one sample per symbol.
% The pulses can be normalised towards a unit peak value or a unit energy
% (in the absence of modulation).
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_elecmod.pulse_shape = 'rc';%'rrc';'sinc';'sech';'gaussian';
% params_elecmod.roll_off = 1;
% params_elecmod.symbol_rate = symbol_rate;
% params_elecmod.fwhm = 1/symbol_rate/6;
% params_elecmod.order = 1;
% params_elecmod.normalisation = 'peak';
% sig = elec_modulator(symbols,params_elecmod); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% symbs             complex symbols to be modulated [complex vector]
%
% params            modulation parameters [structure]
%
%                       params.pulse_shape
%                           shape of the pulses [string]
%
%                           params.pulse_shape ='rc': raised-cosine
%                           params.pulse_shape ='rrc': root raised-cosine
%                           params.pulse_shape ='sinc': sinc 
%                           params.pulse_shape ='sech: sech
%                           params.pulse_shape ='gauss': Gaussian
%
%                       params.symbol_rate
%                           symbol rate of the signal [real scalar]
%
%                           Used for 'rc', 'rrc' and 'sinc' pulses.
%
%                       params.roll_off
%                           roll-off factor [real scalar]
%
%                           Takes values within [0,1]
%                           Used for 'rc' and 'rrc' pulses.
% 
%                       params.fwhm
%                           pulse full-width at half-maximum [real scalar]
%
%                           Used for 'sech' and 'gaussian' pulses.
%
%                       params.order
%                           pulse order [integer scalar]
%
%                           Uses for 'gaussian' pulses.
% 
%                       params.normalisation
%                           pulse normalisation [string]
%
%                           params.normalisation = 'peak'
%                               The pulses are normalised so that their
%                               peak amplitude is one (in the absence of
%                               modulation)
%                           params.normalisation = 'energy'
%                               The pulses are normalised so that their
%                               energy is equal to 1 (in the absence of
%                               modulation)c
%                           params.normalisation = 'default'
%                               In some cases we want to live with the
%                               defaut normalisation of the pulse shape, as
%                               defined in the pulse shape function.
%                               For instance the raised-cosine pulse shape
%                               has unit peak while the root-raised-cosine
%                               pulse shape has unit energy (as defined in
%                               their respective functions; and there are 
%                               good reasons to be so). 
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               modulated electrical signal [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% time_array            time samples, in s [real vector]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


global time_array 


nsamples = length(time_array);
% Number of time samples
nsymbols = length(symbs);
% Number of symbols

if rem(nsamples,nsymbols) == 0
    % Check that the number of samples is an integer multiple of the number
    % of symbols
    samples_per_symbol = nsamples/nsymbols;
    % Number of samples per symbol
else
    error('elec_modulator: the number of samples should be an integer multiple of the number of symbols.');
end


t0 = time_array(length(time_array)/2);
% Position of the isolated pulses that will be used in the convolution with
% the Dirac comb.


% First we generate the isolated pulses shapes.
switch params.pulse_shape
    
    case 'rc'
        % Raised-cosine pulse
        
        ts = 1/params.symbol_rate;        
        sig = elec_pulse_rc(time_array,t0,ts,params.roll_off);       
        
        switch params.normalisation            
            case 'peak'                
                signorm = 1;                
            case 'energy'                
                signorm = ts*(4 - params.roll_off)/4;    
            case 'default'
                signorm = 1;
        end
        
    case 'rrc'
        % Root-raised-cosine pulse
        
        ts = 1/params.symbol_rate;        
        sig = elec_pulse_rrc(time_array,t0,ts,params.roll_off);
        
        switch params.normalisation            
            case 'peak'                
                signorm = (1 - params.roll_off + 4*params.roll_off/pi)/sqrt(ts);                
            case 'energy'                
                signorm = 1;          
            case 'default'
                signorm = 1;
        end
        
    case 'sinc'
        % sinc pulse
        
        ts = 1/params.symbol_rate;
        
        sig = elec_pulse_sinc(time_array,t0,ts);
        
        switch params.normalisation            
            case 'peak'                
                signorm = 1;                
            case 'energy'                
                signorm = params.ts;   
            case 'default'
                signorm = 1;
        end
        
    case 'sech'
        % sech pulse
        
        sig = elec_pulse_sech(time_array,t0,params.fwhm);
        
        switch params.normalisation            
            case 'peak'                
                signorm = 1;                
            case 'energy'                  
                signorm = 1;   %TO DO   
            case 'default'
                signorm = 1;
        end
        
    case 'gauss'
        % Gaussian pulse
        
        sig = elec_pulse_gauss(time_array,t0,params.fwhm,params.order);      
        
        switch params.normalisation            
            case 'peak'                
                signorm = 1;                
            case 'energy'                  
                signorm = 1;  %TO DO   
            case 'default'
                signorm = 1;    
        end        
        
        
    otherwise
        
        error('elec_modulator: pulse shape not defined.');
        
end

sig = sig/signorm;
% The signal is normalised to a peak amplitude of 1

seq = zeros(1,nsamples);
% Create an empty array for the upsampled sequence
seq(samples_per_symbol/2:samples_per_symbol:end) = symbs;
% Upsample the input data symbols with 1 sample per symbol and zeros
% in-between

sig = ifft(fft(seq).*fft(sig));
% Modulated signal. Convolution of the Dirac comb and the pulse shape

sig = circshift(sig,[0 -nsamples/2+1]);
% Compensate delay

end