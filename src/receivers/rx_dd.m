function sig = rx_dd(sig,params)
% Receiver incl. optical filtering and direct and interferometric detection
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements a general receiver using either direct detection
% (DD) or interferometric detection (ID) with either single-ended or
% balanced detection. The main purpose of this receiver is to be called
% from generic Monte-Carlo or Karhunen-Loeve functions that are independent
% of the modulation format employed.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_rx.type = 'dd';%'id';
% params_rx.obpf.type = 'gaussian';
% params_rx.obpf.order = 1;
% params_rx.obpf.bandwidth = 4*symbol_rate;
% params_rx.obpf.centre_frequency = 0;
%
% % Parameters for DD receiver:
% params_rx.elpf.type = 'bessel';
% params_rx.elpf.order = 4;
% params_rx.elpf.f3dB = 0.75*symbol_rate;
% params_rx.pd.responsivity = 1;
% params_rx.pd.shot_noise = 0;
% params_rx.pd.thermal_noise_density = 0;
% params_rx.pd.dark_current = 0;
%
% % Parameters for ID receiver:
% params_rx.mzdi.input_port = 'upper';%'lower';
% params_rx.mzdi.mode = 'tuned';%'general';
% params_rx.mzdi.delay = 1/symbol_rate;
% params_rx.mzdi.phase_shift = 0;
% 
% params_rx.upper.pd.responsivity = 1;
% params_rx.upper.pd.shot_noise = 0;
% params_rx.upper.pd.thermal_noise_density = 0;
% params_rx.upper.pd.dark_current = 0;
% params_rx.lower.pd.responsivity = 1;
% params_rx.lower.pd.shot_noise = 0;
% params_rx.lower.pd.thermal_noise_density = 0;
% params_rx.lower.pd.dark_current = 0;
% 
% params_rx.upper.elpf.type = 'bessel';
% params_rx.upper.elpf.order = 4;
% params_rx.upper.elpf.f3dB = 0.75*symbol_rate;
% params_rx.lower.elpf.type = 'bessel';
% params_rx.lower.elpf.order = 4;
% params_rx.lower.elpf.f3dB = 0.75*symbol_rate;
% params_rx.bd.detection_mode = 'balanced';%'single_ended_upper';%'single_ended_lower';
% params_rx.bd.electrical_delay = 0;
% params_rx.bd.polarity = 1;%-1.
% 
% sig = rx_dd(sig,params_rx);
% 
% % For MZDI input_port = 'upper' and mode = 'Tuned'
% % 'single_ended_upper' -> compare to original data sequence and ignore 
% % 1st bit. The optical signal is AMI.
% % 'single_ended_lower' -> compare to inverted original data sequence and 
% % ignore 1st bit. The optical signal is DB.
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
%                       
% params            structure containing the receiver parameters
%                       [structure]
%
%                       params.type
%                           type of receiver [string]
%
%                           params.type = 'dd': direct detection receiver
%                           params.type = 'id': interferometric detection 
%                               receiver  
%
%                       params.obpf 
%                           optical bandpass filter parameters [structure]
%
%                           See opt_opbf function for more details.
%
%                           params.obpf.type 
%                           params.obpf.order
%                           params.obpf.bandwidth
%                           params.obpf.centre_frequency
% 
%                       params.pd 
%                           photodiode parameters, for DD receiver 
%                           [stucture]
%
%                           See rx_pin for more details.
%
%                           params.pd.responsivity
%                           params.pd.thermal_noise_density
%                           params.pd.shot_noise
%                           params.pd.dark_current
%
%                       params.elpf
%                           post-detection filter parameters, for DD
%                           receiver [structure]
%
%                           See elec_elpf for more details.
%
%                           params.elpf.type
%                           params.elpf.order
%                           params.elpf.f3dB
% 
%                       params.mzdi
%                           Mach-Zehnder delay interferometer parameters
%                           for ID receiver [structure]
%
%                           See opt_mzdi for more details.
%
%                           params.mzdi.input_port ['upper'/'lower']
%                           params.mzdi.mode ['tuned'/'general']
%                           params_mzdi.coupling_ratio_in
%                           params_mzdi.coupling_ratio_in
%                           params.mzdi.delay_target
%                           params.mzdi.phase_shift
%
%                       params.upper.pd 
%                           upper photodiode parameters, for ID receiver 
%                           [stucture]
%
%                           params.upper.pd.responsivity
%                           params.upper.pd.thermal_noise_density
%                           params.upper.pd.shot_noise
%                           params.upper.pd.dark_current
% 
%                       params.lower.pd 
%                           lower photodiode parameters, for ID receiver 
%                           [stucture]
%
%                           params.lower.pd.responsivity
%                           params.lower.pd.thermal_noise_density
%                           params.lower.pd.shot_noise
%                           params.lower.pd.dark_current
%
%                       params.upper.elpf
%                           upper post-detection filter parameters, for ID
%                           receiver [structure]
%
%                           params.upper.elpf.type
%                           params.upper.elpf.order
%                           params.upper.elpf.f3dB
%
%                       params.lower.elpf
%                           lower post-detection filter parameters, for ID
%                           receiver [structure]
%
%                           params.lower.elpf.type
%                           params.lower.elpf.order
%                           params.lower.elpf.f3dB
%
%                       params.bd
%                           balanced detection parameters, for ID receiver
%                           [structure]
%
%                           params.bd.detection_mode 
%                               Determines whether balanced or single-ended
%                               detection is used [string]
%
%                               params.bd.detection_mode = 'balanced';
%                               params.bd.detection_mode = 'single_ended_upper';
%                               params.bd.detection_mode = 'single_ended_lower';
%
%                           params.bd.electrical_delay
%                               Delay between the upper and lower arm 
%                               photocurrents before they are subtracted, 
%                               in s [real scalar]
%
%                           params.bd.polarity
%                               Polarity of the balanced-detected signal
%                               [1,-1]
%
%                               params.bd.polarity = 1; (normal)
%                               params.bd.polarity = -1; (inverted)
%                               Possibility to invert the balanced detected 
%                               signal. Will be required for 67% RZ DPSK. 
%                               Otherwise use normal for all the other DPSK 
%                               formats.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               optically filtered, detected and electrically filtered 
%                       signal [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% frequency_array       relative frequency samples, in Hz [real vector]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global frequency_array

otf = opt_tf_obpf(params.obpf,frequency_array);
sig = opt_filter(sig,otf);
% Optical bandpass filtering of the received optical signal

if strcmp(params.type,'dd')
    % Direct detection receiver.
    
    sig = rx_pin(sig,params);
    % Photocurrent generated by the photodiode, including low-pass
    % filtering

elseif strcmp(params.type,'id')
    % Interferometric detection receiver
    
    if strcmp(params.mzdi.input_port,'upper')
        sig11 = sig;
        sig12 = opt_nosig;

    elseif strcmp(params.mzdi.input_port,'lower')
        sig11 = opt_nosig;
        sig12 = sig;

    else
        disp('rx_dd: non existing MZDI input port.');
    end
    % Select input port of MZDI    


    [sig21,sig22] = opt_mzdi(sig11,sig12,params.mzdi);
    % Mach-Zehnder delay interferometer

    i_upper = rx_pin(sig21,params.upper);
    % Detection and low-pass filtering, upper photodiode.
    i_lower = rx_pin(sig22,params.lower);
    % Detection and low-pass filtering, lower photodiode
    
    if strcmp(params.bd.detection_mode,'balanced')
        % Balanced detection
        
        [ilowerdelayed,~] = elec_delay(i_lower,params.bd.electrical_delay);
        sig = i_upper - ilowerdelayed;
        
    elseif strcmp(params.bd.detection_mode,'single_ended_upper')
        % Single-ended detection, upper arm.
        
        sig = i_upper;
        
    elseif strcmp(params.bd.detection_mode,'single_ended_lower')
        % Single ended detection, lower arm.
        
        sig = i_lower;
        
    else
        disp('rx_dd: balanced detection mode not implemented.');
    end
    % Return photocurrent, depending on whether balanced or single ended
    % detection is used.  
    
    sig = sign(params.bd.polarity)*sig;
    % Inversion of the signal, if necessary.
    
else
    disp('rx_dd: receiver type not implemented.');
end

end