function sig = opt_laser_cw(params)
% Ideal continuous wave laser
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates an optical continuous wave (CW) signal with a 
% specified emission frequency, power and linewidth. 
% The signal is polarised along the -x direction.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_cw.power = 1.0e-3;
% params_cw.linewidth = 0;
% params_cw.emission_frequency = 193.1e12;
% sig = opt_laser_cw(params_cw);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% params            CW laser parameters [structure]
%   
%                       params.power 
%                           signal average power, in W [real scalar]
%
%                       params.linewidth
%                           signal linewidth, in Hz [real scalar]
%
%                       params.emission_frequency
%                           desired laser emission frequency, in Hz
%                           [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output optical signal [optical signal structure]
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
% df                    frequency samples separation, in Hz [real scalar]
%
% -------------------------------------------------------------------------
% REMARKS:
% -------------------------------------------------------------------------
% The structure params is updated with the fields:
%
% params.actual_emission_frequency      actual laser emission frequency, 
%                                           in Hz [real scalar]
%
% params.frequency_offset               frequency offset between the actual
%                                           emission frequency and the
%                                           desired emission frequency,
%                                           in Hz [real scalar]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global reference_frequency 
global time_array 
global dt
global df 

params.actual_emission_frequency = reference_frequency + round((params.emission_frequency - reference_frequency)/df)*df;
% Ensure that the emission frequency falls onto the frequency grid

params.frequency_offset = params.actual_emission_frequency - params.emission_frequency;
% Frequency offset between the actual emission frequency and the desired
% emission frequency

fprintf(1,'\n\n%s\n%s\t%f\t%s\n','opt_laser_cw:','Desired emission frequency: ',params.emission_frequency/1.0e12,' THz');
fprintf(1,'%s\t\t%f\t%s\n','Actual emission frequency: ',params.actual_emission_frequency/1.0e12,' THz');
fprintf(1,'%s\t\t\t\t\t\t%f\t%s\n','Delta f: ',params.frequency_offset/1.0e9,' GHz');
% Notify the user that the emission frequency has been adjusted

nsamples = length(time_array);
% Number of samples in the signal to generate

wiener_phase_noise = cumsum(sqrt(2*pi*params.linewidth*dt)*randn(1,nsamples));
% Wiener phase noise

wiener_phase_noise = wiener_phase_noise - wiener_phase_noise(1);
% We remove a constant phase shift so that the first value of the phase
% noise is zero.


% Apply correction factor to the phase noise process in order to avoid
% phase discontinuities in conjunction with the periodic nature of the
% signal.
% C. J. Rasmussen, "Transmission analysis in WDM networks", Ph.D. Thesis,
% Technical University of Denmark, 1999.
% pp. 118-119.
phase_diff = wiener_phase_noise(nsamples) - wiener_phase_noise(1);
if abs(phase_diff) > eps
    % No need to worry if the linewidth is zero or if the first and last
    % samples have the same phase (the phase may be continuous though with
    % phase_diff ~= 0 since phase(1) = phase(nsamples +1) in this case.
    delta_phi = round(phase_diff/(2*pi))*2*pi;
    phase_correction = time_array*(delta_phi - phase_diff)/((nsamples - 1)*dt);
    wiener_phase_noise = wiener_phase_noise + phase_correction;
end


sig = struct;
sig.x = sqrt(params.power).*exp(-1i*wiener_phase_noise).*exp(1i*2*pi*(params.actual_emission_frequency - reference_frequency)*time_array);
sig.y = zeros(1,nsamples);
% Returns the complex field

end