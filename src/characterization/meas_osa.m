function varargout = meas_osa(sig,params)
% Optical spectrum analyser
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This functions visualises optical spectra averaged over a given bandwidth 
% and save those into files.
% This is not a real implementation of an optical spectrum analyser. The
% function calculates the power spectrum of the input signal and performs
% some smoothing via adjacent averaging.
% The resolution_bandwidth parameter is actually rather a video bandwidth,
% even though the desired smoothing effect will take place with values
% similar to those normally used for the resolution bandwidth.
% This implementation is sufficient for most practical uses.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_osa.pol = 'x';%'y';'both';
% params_osa.display_interval = [frequency_array(1) frequency_array(end)];
% params_osa.resolution_bandwidth = 0;%12.5e9;
% params_osa.sensitivity = -40;
% params_osa.display = 1;
% params_osa.save.txt = 0;
% params_osa.save.emf = 0;
% params_osa.save.jpg = 0;
% params_osa.name = 'Optical spectrum';
% meas_osa(sig.x,params_osa); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               optical signal to be characterised 
%                       [optical signal structure]
%
% params            OSA parameters [structure]
%
%                       params.pol
%                           polarisation of the signal that will be 
%                           visualised [string]
%
%                           params.pol = 'x'
%                           params.pol = 'y'
%                           params.pol = 'both'
%
%                       params.display_interval
%                           frequency range over which the spectrum will be
%                           displayed, in Hz [2-element vector]
%
%                           params.display_interval = [start stop]
%
%                       params.resolution_bandwidth
%                           resolution bandwidth, in Hz [real scalar]
% 
%                           Adjacent averaging filtering is performed on
%                           the power spectrum, therefore this parameter is
%                           technically rather a video bandwidth.
%
%                           Rectangular filtering is applied. 
%
%                           If params.resolution_bandwidth = 0, no 
%                           resolution bandwidth is applied.
%
%                       params.sensitivity
%                           sensitivity of the optical spectrum analyser, 
%                           in dBm [real scalar]
%
%                           Power values below the sensitivity level will 
%                           not be displayed.
%
%                       params.display
%                           toggles the display on or off [0/1]
%
%                           params.display = 0;
%                           params.display = 1;
%
%                       params.name
%                           name of the visualiser and the saved files
%                           [string]
%
%                       params.save = 
%                           saving information [structure]
%
%                               params.save.txt
%                                   switch determining whether to save the 
%                                   spectrum in a .dat text file [0/1] 
%
%                               params.save.emf
%                                   switch determining whether to save the 
%                                   trace of the spectrum in an .emf file,
%                                   without legend nor axes, in order to be
%                                   directly used in e.g. pptx 
%                                   presentations) [0/1] 
%
%                               params.save.jpg
%                                   switch determining whether to save a 
%                                   picture of the displayed spectrum in a 
%                                   .jpg file [0/1]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% varargout         optional export of spectrum data [structure]
%
%                       varargout(1)
%                           spectrum in linear scale, in W [real vector]
%
%                       varargout(2)      
%                           smoothed spectrum in linear scale, in W 
%                           [real vector]
%
%                           Note that the effect of the smoothing filter at
%                           the edge of the frequency range is not 
%                           corrected for. 
%                           The length varargout(2) is the same as the 
%                           length of frequency_array. 
%                           One may need to truncate the spectrum for nice
%                           looking display, as is done in the part of 
%                           this function dealing with the display of the 
%                           spectrum.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% frequency_array   relative frequency samples, in Hz [real vector]
%
% df                frequency samples separation, in Hz [real scalar]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global df 
global frequency_array


params.display_interval = params.display_interval(:).';
% Force display_interval to a line vector

if params.display_interval(1) > params.display_interval(2)
    fliplr(params.display_interval);
end
% Ensure the display interval is in increasing value order

nsamples = length(sig.x);
% Number of samples in the input signal

if strcmp(params.pol,'x')
    sig_test = sig.x;
elseif strcmp(params.pol,'y')
    sig_test = sig.y;
elseif strcmp(params.pol,'both')
    sig_test = sqrt(abs(sig.x).^2 + abs(sig.y).^2).*exp(-1j*0.5*(angle(sig.x) + angle(sig.y)));
end
% Signal to plot, depending on the polarisation selection

sig_spectrum = fftshift(fft(sig_test))/nsamples;
% Calculate the spectrum from the time domain signal

sig_power = abs(sig_spectrum).^2;
% Power spectrum

varargout(1) = {sig_power};
% Save the spectrum in linear scale as first optional output argument

sig_power_log = 10*log10((sig_power + 1e-200)/1e-3);
% Power in dBm


% In the following block of code, averaging over the given resolution
% bandwidth is performed.
if params.resolution_bandwidth == 0
    sig_smoothed = sig_power_log;
    shift = 0;   
    
    varargout{2} = sig_power;
    
else
    samples = round(params.resolution_bandwidth/df);
    % Calculates the number of samples corresponding to a given resolution
    % bandwidth
    if rem(samples,2) == 1
        samples = samples - 1;
    end
    % Adjust this number to an even number in case round(params.resolution_bandwidth/df) 
    % turns out to be odd. Required for proper compensation of the shift 
    % induced by the use of the filter function.

    sig_smoothed = filter(ones(1,samples),1,sig_power);
    shift = -samples/2;
    % Perform sum of the values of the power samples falling within the
    % resolution bandwidth. Calculate values of the shift to apply back to 
    % the signal due to the use of the filter function.

    % smoothed=filter(ones(1,samples)/samples,1,power);
    % shift=-samples/2;
    % this code if adjacent avergaging only (no power integration in the
    % resolution bandwidth) is performed

    sig_smoothed = circshift(sig_smoothed',shift)';
    % Compensate for the sample shift induced by the use of the filter
    % function
    
    varargout(2) = {sig_smoothed};
    % Save the smoothed spectrum in linear scale as optional output
    % argument
    
    sig_smoothed = 10*log10(sig_smoothed/1e-3);
    % Convert the smoothed spectrum into log scale (dBm)
    
end
% End of integration of the spectrum over a specified resolution bandwidth.


frequency_array_local = frequency_array(1:nsamples + shift);
sig_power_log = sig_power_log(1:nsamples + shift);
sig_smoothed = sig_smoothed(1:nsamples + shift);
% Restrict the frequency range to avoid visualising and saving power that
% has leaked from the low frequency to the high frequency region due to the
% use of circshift to compensate for frequency offset.

if isempty(params.display_interval)
    freq_min = frequency_array_local(1);
    freq_max = frequency_array_local(length(frequency_array_local));
else
    freq_min = params.display_interval(1);
    freq_max = params.display_interval(2);
end
power_min = params.sensitivity;
power_max = max(sig_smoothed) + 10;

if power_min >= power_max
    disp('meas_osa: sensitivity of the OSA is too high. It will be decreased by 20 dB.' );
    power_min = power_max - 20;
end

if params.display == 1
    % Display spectrum only when asked explicitely
    
    optical_spectrum = figure('Name',params.name);
    plot(frequency_array_local./1.0e9,sig_smoothed,'Color',[0 0 1],'LineWidth',1);
    axis([freq_min/1.0e9 freq_max/1.0e9 power_min power_max]);
    grid on;
    xlabel('relative frequency (GHz)');
    ylabel('power (dBm)');
    
end

savearray = [frequency_array_local',sig_power_log',sig_smoothed'];
% Build array containing the information to be saved in the txt file


if params.save.jpg == 1
    file_name_jpg = [params.name,'.jpg'];    
    print(optical_spectrum,'-djpeg','-r300',file_name_jpg); 
end
% Save into jpg file

if params.save.txt == 1
    file_name_txt = [params.name,'.dat'];    
    fid = fopen(file_name_txt,'w');
    fprintf(fid,'%s\t%s\t%s\n','freq','spectrum','smoothed');
    fprintf(fid,'%s\t%s\t%s\n','(Hz)','(dBm)','(dBm)');
    fprintf(fid,'%g\t%g\t%g\n',savearray');
    fclose(fid);
end
% Save into txt file


if params.save.emf == 1
    file_name_emf = [params.name,'.emf'];
    temp = figure();
    plot(frequency_array_local./1.0e9,sig_smoothed,'LineWidth',1,'Color',[0 0 0]);
    axis([freq_min/1.0e9 freq_max/1.0e9 power_min power_max]);
    axis off;
    print(temp,'-dmeta','-r300',file_name_emf);
    close(temp);
end
% Save the integrated spectrum into emf file

end