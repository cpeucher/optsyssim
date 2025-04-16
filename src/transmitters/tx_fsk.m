function sig = tx_fsk(bit_pattern,params)
% Ideal optical binary FSK transmitter
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements an ideal optical binary frequency shift keying 
% (FSK) transmitter enabling the generation of FSK (non-continuous phase),
% continuous-phase frequency shift keying (CPFSK), minimum shift keying 
% (MSK) and Gaussian minimum shift keying (GMSK) formats.
% The FSK signal is generated from two sinusoidal signals at different
% frequencies, while the continuous phase formats (CPFSK, MSK, GMSK)
% are obtained through the generation of a continuous phase term by
% integration of the modulating signal. In case of MSK, the modulation
% index, defined as h=\Delta f / Rs, where \Delta f is the peak to peak
% frequency deviation and Rs is the bit rate is equal to 0.5. For GSMK
% modulation, the NRZ electrical signal is filtered by a Gaussian filter
% with specified normalised bandwidth BT before being used for MSK
% modulation.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_tx.type = 'fsk';%'cpfsk';'msk';'gmsk';
% params_tx.symbol_rate = symbol_rate;
% params_tx.emission_frequency = reference_frequency;
% params_tx.tone_spacing = 2*params_tx.symbol_rate; % For FSK and CPFSK only.
% params_tx.rise_time = 1/params_tx.symbol_rate/4; 
% params_tx.bt = 0.3;                               % For GMSK only.
% params_tx.power = 1.0e-3;
% sig = tx_fsk(bit_pattern,params_tx);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% bit_pattern       binary bit pattern to be modulated [binary vector]
%
% params            FSK transmitter parameters [structure]
%
%                       params.type
%                           type of modulation [string]
%
%                           params.type = 'fsk';
%                           params.type = 'cpfsk';
%                           params.type = 'msk';
%                           params.type = 'gmsk';
%
%                       params.emission_frequency
%                           transmitter emission frequency, in Hz
%                           [real scalar]
%
%                           This corresponds to the centre frequency of the
%                           spectrum, i.e. the modulation frequency tones
%                           are located symetrically around this value.
%
%                       params.tone_spacing
%                           FSK tone spacing, in Hz [real scalar]
%
%                           The modulation index is defined as 
%                           h = params.tone_spacing/params.symbol_rate.
%
%                           This value is only required for FSK and CPFSK. 
%                           For MSK and GMSK it is automatically adjusted
%                           to h = 0.5, i.e. 
%                           params.tone_spacing = 0.5*params.symbol_rate.  
%
%                       params.symbol_rate
%                           symbol rate of the generated signal, in bit/s
%                           [real scalar]
%
%                       params.rise_time
%                           rise time of the electrical signal generating 
%                           the frequency modulation, in s [real scalar]
%
%                           This value is required for FSK, CPFSK, MSK.
%                           In GMSK, this parameter is determined from 
%                           params.bt.
%
%                       params.bt
%                           normalised bandwidth of the Gaussian filter
%                           in the GMSK case [real scalar]
%
%                           The bandwidth is normalised to
%                           the symbol rate.
%
%                       params.power
%                           signal power, in W [real scalar]
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
% dt                    time samples separation, in s [real scalar]
%
% reference_frequency   reference frequency, in Hz [real scalar]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global reference_frequency
global time_array
global dt


switch params.type
    
    case 'msk'
        % MSK modulation
        
        data_sig = elec_pulse_sequence_nrz(bit_pattern,params.rise_time);
        % Generation of a rectangular pulse sequence with rise time        
  
        data_sig = (data_sig - min(data_sig))/(max(data_sig) - min(data_sig));
        % It is ensured the modulating signal is within [0 1]
        
        data_sig = (2*data_sig - 1)/2;
        % The data signal is now bipolar and within [-1/2 1/2]
        
        params.tone_spacing = params.symbol_rate/2;
        % For MSK the tone spacing is equal to half the symbol rate
        
        sig_phase = 2*pi*params.tone_spacing*cumtrapz(data_sig)*dt;
        % Integrate the data waveform to generate the continuous phase
        
        
    case 'gmsk'   
        % GMSK modulation
        
        data_sig = elec_pulse_sequence_nrz(bit_pattern,0);
        % Generation of a rectangular pulse sequence
        
        rise_time_filter.type = 'gaussian';
        rise_time_filter.order = 1;
        rise_time_filter.f3dB = params.bt*params.symbol_rate;
        % Filtering by a Gaussian ELPF whose normalised f3dB is params.bt        
        
        data_sig = elec_elpf(data_sig,rise_time_filter);
        % A Gaussian low-pass filter is used to introduce the rise time
        
        data_sig = (data_sig - min(data_sig))/(max(data_sig) - min(data_sig));
        % It is ensured the modulating signal is within [0 1]
        
        data_sig = (2*data_sig - 1)/2;
        % The data signal is now bipolar and within [-1/2 1/2]
                
        params.tone_spacing = params.symbol_rate/2;
        % For MSK the tone spacing is equal to half the symbol rate
        
        sig_phase = 2*pi*params.tone_spacing*cumtrapz(data_sig)*dt;
        % Integrate the data waveform to generate the continuous phase
        
        
    case 'cpfsk'
        % CPFSK modulation
        
        data_sig = elec_pulse_sequence_nrz(bit_pattern,params.rise_time);
        % Generation of a rectangular pulse sequence with rise time
        
        data_sig = (data_sig - min(data_sig))/(max(data_sig) - min(data_sig));
        % It is ensured the modulating signal is within [0 1]
        
        data_sig = (2*data_sig - 1)/2;
        % The data signal is now bipolar and within [-1/2 1/2] 
        
        sig_phase = 2*pi*params.tone_spacing*cumtrapz(data_sig)*dt;
        % Integrate the data waveform to generate the continuous phase
        
    case 'fsk'
        % FSK modulation.
        
        data_sig = elec_pulse_sequence_nrz(bit_pattern,params.rise_time);
        % Generation of a rectangular pulse sequence with rise time
        
        data_sig = (data_sig - min(data_sig))/(max(data_sig) - min(data_sig));
        % It is ensured the modulating signal is within [0 1]
        
        data_sig = (2*data_sig - 1)/2;
        % The data signal is now bipolar and within [-1/2 1/2]
        
        sig_phase = 2*pi*params.tone_spacing*data_sig.*time_array;
        % Integrate the data waveform to generate the continuous phase 
        
        
    otherwise
        
        error('tx_fsk: modulation type not implemented.');

end

sig.x = sqrt(params.power)*exp(1i*2*pi*(params.emission_frequency - reference_frequency).*time_array + 1i*sig_phase);
sig.y = zeros(1,length(time_array));
% Optical signal

end