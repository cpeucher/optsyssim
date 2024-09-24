function  [sig,line_indices] = opt_comb(params)
% Arbitrary frequency comb generation
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates an arbitray frequency comb at specified
% frequencies. The power, polarisation, initial phase or linewidth of each
% line in the comb can be specified individually. The lines may be
% completely coherent, or they may be impacted by individual realisations
% of phase noise. The comb may be used as e.g. broadband source, test
% signal for amplifier characterisation, WDM source etc.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_comb.nlines = 13;
% % Number of lines in the comb.
% params_comb.centre_frequency = reference_frequency;
% % Centre frequency of the comb, in Hz.
% params_comb.line_spacing = 100e9;
% % Frequency spacing of the lines in the comb, in Hz.
% params_comb.frequency = params_comb.centre_frequency - (params_comb.nlines - 1)/2*params_comb.line_spacing+(0:1:params_comb.nlines - 1)*params_comb.line_spacing;
% % Frequencies of the lines in the comb, in Hz.
% params_comb.power = 1e-3*randn(1,params_comb.nlines);
% % Power of the lines in the comb, in W.
% params_comb.phase = zeros(1,params_comb.nlines);
% % Phases of the lines in the comb, in rad.
% params_comb.linewidth = 1e6*rand(1,params_comb.nlines);
% % Linewidths of the lines in the comb, in Hz.
% params_comb.jones.x = ones(1,params_comb.nlines);
% params_comb.jones.y = zeros(1,params_comb.nlines);
% % Jones vectors of the lines in the comb. Should be normalised.
% [sig,params_comb.line_indices] = opt_comb(params_comb);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% params            frequency comb parameters [structure]
%
%                       params.power
%                           power of each line in the comb, in W
%                           [real vector]
%
%                       params.frequency
%                           frequencies of the comb lines, in Hz
%                            [real vector]
%
%                       params.phase
%                           phases of the comb lines, in radians
%                           [real vector]
%
%                       params.linewidth
%                           linewidths the comb lines, in Hz
%                           [real vector]
%
%                           Obviously the params.phase vector is not very 
%                           meaningful when the linewidths are non zero.
%                           The initial phases params.phase are to be used 
%                           primarily to set the relative phases of a comb 
%                           with infinitely narrow lines.
%
%                       params.jones
%                           definition of the state of polarisation of 
%                           each line in the comb [structure]
%
%                               params.jones.x
%                                   -x component of the normalised Jones 
%                                   vector [complex vector]
%
%                               params.jones.y
%                                   -y component of the normalised Jones 
%                                   vector [complex vector]
%
%                                   The Jones vector should be normalised,
%                                   i.e. abs(jones.x)^2 + abs(jones.y)^2 =1
%                                   for each of the lines.
%
%                       Note that all input arrays should have the same 
%                       length.
%                       The frequencies of the comb can be arbitrary, i.e. 
%                       not necessarily regularly spaced. 
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output optical comb signal [optical signal structure]
%
% line_indices      indices of frequency_array corresponding to the
%                       lines in the comb [integer vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% time_array            time samples, in s [real vector]
%
% dt                    time samples separation, in s [real scalar]
%
% frequency_array       relative frequency samples, in Hz [real vector]
%
% reference_frequency   reference frequency, in Hz [real scalar]
%
% df                    frequency samples separation, in Hz [real scalar]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global time_array 
global frequency_array 
global dt 
global df 
global reference_frequency



nsamples = length(time_array);
% Number of samples in the signal

try
    input_test = params.power.*params.frequency.*params.phase.*params.linewidth.*params.jones.x.*params.jones.y;
catch
    error('opt_comb: input arrays should all have the same length.');
end
% Check that all input arrays have the same length

sig.x = zeros(1,nsamples);
sig.y = zeros(1,nsamples);
% Initialise the output signal


for iline = 1:length(params.frequency)
    % Loop over the number of lines in the comb.
    
    params.frequency(iline) = reference_frequency + round((params.frequency(iline)-reference_frequency)/df)*df;
    % Ensure the emission frequency falls onto the frequency grid  
    
    line_indices(iline) = find(frequency_array + reference_frequency == params.frequency(iline));
    % Locate the indices of the centre frequency of the line
    
    field = sqrt(params.power(iline)).*exp(-1i*params.phase(iline)).*exp(1i*2*pi*(params.frequency(iline) - reference_frequency)*time_array);
    % Power, phase and frequency terms of the complex field
    
    % We generate the phase noise realisation corresponding to
    % the linewidth specified for this given spectral line. 
    if params.linewidth(iline) == 0
        phase_noise = zeros(1,nsamples);
        % If the lines are infinitively narrow, the phase noise realisation
        % reduces to 0
    else
        phase_noise = cumsum(sqrt(2*pi*params.linewidth(iline)*dt)*randn(1,nsamples));
        % Wiener phase noise
        
        % Is what follows necessary in this case ? 
        % Apply correction factor to the phase noise process in order to avoid
        % phase discontinuities in conjunction with the periodic nature of the
        % signal.
        % C. J. Rasmussen, "Transmission analysis in WDM networks", Ph.D. thesis,
        % Technical University of Denmark, 1999.
        % pp. 118-119.
        phase_diff = phase_noise(nsamples) - phase_noise(1);
        if phase_diff ~= 0
            % No need to worry if the linewidth is zero or if the first and last
            % samples have the same phase (the phase may be continuous though with
            % PhaseDiff ~= 0 since phase(1)=phase(Nsamples +1) in this case.
            delta_phi = round(phase_diff/(2*pi))*2*pi;
            phase_correction = time_array*(delta_phi - phase_diff)/((nsamples - 1)*dt);
            phase_noise = phase_noise + phase_correction;
        end
    end
    % End of phase noise realisation generation.    
    
    sig.x = sig.x + field.*exp(-1i*phase_noise)*params.jones.x(iline);
    sig.y = sig.y + field.*exp(-1i*phase_noise)*params.jones.y(iline);
    % Add the new line to the comb  

end
% End of loop over lines in the frequency comb

end