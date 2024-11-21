function spectrum = meas_esa(sig,params)
% Electrical spectrum analyser
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function displays the spectra of electrical signals and therefore
% emulates an electrical spectrum analyser (ESA).
% The single-sided spectrum is represented, i.e. the total power can be
% obtained by integration over the positive half-frequency axis.
% An ideal rectangular filter is used to integrate the power within the
% resolution bandwidth.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_esa.display_interval = [0 frequency_array(end)];
% params_esa.resolution_bandwidth = 0;
% params_esa.input_impedance = 1;
% params_esa.display = 1;
% params_esa.save.txt = 0;
% params_esa.save.emf = 0;
% params_esa.name = 'RF spectrum';
% spectrum = meas_esa(sig,params_esa); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               electrical signal to be characterised [real vector]
%
% params            radio frequency spectrum display parameters
%                       [structure]
%
%                       params.display_interval
%                           frequency interval over which the RF spectrum 
%                           will be displayed, in Hz 
%                           [2-element real vector]
%
%                           params.display_interval = [start stop];
%
%                       params.resolution_bandwidth
%                           resolution bandwidth of the displayed RF 
%                           spectrum, in Hz [real scalar]
%
%                       params.input_impedance
%                           input resistance of the ESA, in ohms
%                           [real scalar]
%
%                       params.display
%                           specifies whether the RF spectrum is displayed 
%                           in a figure, or just calculated and 
%                           returned in the relevant output structure [0/1]
%
%                           params.display = 1;
%                           params.display = 0;
%
%                       params.name
%                           name of the visualiser and output files
%                           [string]
%
%                       params.save 
%                           saving information [structure]
%
%                           params.save.txt
%                               switch to to save spectrum in .dat text
%                               file [0/1]
%
%                           params.save.emf
%                               switch to save the trace of the spectrum in
%                               emf file  (without legend nor axes, in 
%                               order to be directly used in e.g. pptx 
%                               presentation) [0/1]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% spectrum          calculated RF spectrum [structure]
%
%                       Only the positive frequencies over the entire 
%                       simulation bandwidth are retained.
%                       The spectra are one-sided.
%
%                       spectrum.frequency
%                           positive RF frequencies over which the raw and 
%                           smoothed spectra are saved, in Hz
%                           [real vector]
%
%                       spectrum.raw
%                           raw spectrum, in W [real vector]
%
%                       spectrum.smoothed
%                           smoothed (over the resolution bandwidth) 
%                           spectrum, in W [real vector]
%
%                       spectrum.rb
%                           actual resolution bandwidth in which the
%                           smoothed spectrum has been calculated, in Hz
%                           [real scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% df                frequency samples separation, in Hz [real scalar]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global df


% -------------------------------------------------------------------------
% Check input
% -------------------------------------------------------------------------
params.display_interval = params.display_interval(:).';
% Force display_interval to a line vector

if params.display_interval(1) > params.display_interval(2)
    fliplr(params.display_interval);
end
% Ensure the display interval is in increasing value order

if params.display_interval(1) < 0
    error('meas_esa: only positive frequencies are displayed (one-sided power spectrum)');
end
% Ensure the lower visualisation frequency is positive (>= 0)


% -------------------------------------------------------------------------
%  Calculate two-sided spectrum
% -------------------------------------------------------------------------
nsamples = length(sig);
% Number of samples in the electrical input signal

xf2 = fft(sig)/nsamples/sqrt(params.input_impedance);
% fft of the input signal, in fft order
% Also take load resistance into account


% -------------------------------------------------------------------------
% Calculate one-sided spectrum
% -------------------------------------------------------------------------
xf1 = xf2(1:nsamples/2);
% Restrict to positive frequency range (incl. dc) 
xf1(2:end) = sqrt(2)*xf1(2:end);
% Multiply by sqrt(2) apart from the DC term, so that the single-sided 
% spectrum is considered

xf1 = abs(xf1).^2;
% Now power spectrum, in W

freq1 = [0:nsamples/2-1]*df;
% Positive frequencies range

display_indices = find(freq1 >= params.display_interval(1) & freq1 <= params.display_interval(2));
% Find indices corresponding to the display frequency interval


% -------------------------------------------------------------------------
%  Resolution bandwidth
% -------------------------------------------------------------------------
if params.resolution_bandwidth ~= 0
    % If a non-zero value of resolution bandwidth is provided, we need to
    % average the spectrum over this bandwidth. 
    % This needs to be done on the double-sided spectrum.
    
    samples_in_resolution_bandwidth = round(params.resolution_bandwidth/df);
    % Calculate the number of samples corresponding to a given resolution
    % bandwidth
    
    if rem(samples_in_resolution_bandwidth,2) == 1
        samples_in_resolution_bandwidth = samples_in_resolution_bandwidth - 1;
    end
    % Adjust this number to an even number in case round(rb/df) turns out
    % to be odd. Required for proper compensation of the shift induced by
    % the use of the filter function.
    % We should ensure that we still have a sample...
    
    % Pre-process the signal to filter in order to obtain the correct power
    % levels in the averaged one-sided spectrum  
    xf_to_filter = abs(xf2).^2;
    % Power spectrum, in fft order
    xf_to_filter(2:end) = 2*xf_to_filter(2:end);
    % We anticipate the fact that we want the one-sided spectrum at the
    % end, while we perform filtering on the two-sided spectrum.
    % Therefore we double the power of all frequency components apart from 
    % DC. In principe we should not double the power of the frequency
    % component at Nyquist (-fs/2) too. Normally we do not expect the power
    % to be very high there...
    xf_to_filter = fftshift(xf_to_filter);
    % Now in increasing frequency order
    xf_to_filter(1) = xf_to_filter(1)/2;
    % Power at Nyquist frequency should not be doubled
    

    xfs1 = filter(ones(1,samples_in_resolution_bandwidth),1,xf_to_filter);
    % Integrate the power within specified resolution bandwidth
    % Perform sum of the values of the power samples falling within the
    % resolution bandwidth. Calculate values of the shift to apply back to 
    % the signal due to the use of the filter function.
    % xfs1 is in increasing frequency order.
    % xfs1 is power spectrum, in W.
    
    xfs1 = circshift(xfs1',-samples_in_resolution_bandwidth/2)';
    % Compensate for the sample shift induced by the use of the filter
    % function
    % Necessary so that the lines of the spectrum are centered at the
    % correct frequencies.
    
    xfs1 = xfs1(nsamples/2+1:end);
    % Select the positive frequencies only


end
% End of averaging over the specified resolution bandwidth

% We now have:
% xf1       "raw" one-sided spectrum, in W
% xfs1      "averaged" one-sided spectrum, in W
% freq1     corresponding positive frequencies axis


% -------------------------------------------------------------------------
% What to plot and return?
% -------------------------------------------------------------------------
if params.resolution_bandwidth == 0
    % Display raw spectrum
    % Return raw spectrum only
    
    xf_display = xf1;
    
    spectrum.raw = xf1;
    spectrum.smoothed = xf1;
    spectrum.rb = 0;
    
else
    % Display spectrum in specified resolution bandwidth
    % Return raw spectrum + spectrum in specified resolution bandwidth
    
    xf_display = xfs1;
    
    spectrum.raw = xf1;
    spectrum.smoothed = xfs1;
    spectrum.rb = (samples_in_resolution_bandwidth - 1)*df;
    
end
% xf_display is the power spectrum
% Save spectrum data in output structure. Only the positive frequencies
% over the entire simulation bandwidth are retained.

spectrum.frequency = freq1;
% Also return positive frequencies axis


% -------------------------------------------------------------------------
% Display RF spectrum.
% -------------------------------------------------------------------------
if params.display == 1

    figure('Name',params.name);
    plot(freq1(display_indices)/1.0e9,10*log10(xf_display(display_indices)/1.0e-3),'Color',[0 0 1]);
    grid;
    xlabel('frequency (GHz)');
    ylabel('RF power (dBm)');

end
% Display the RF spectrum 
  

% -------------------------------------------------------------------------
% Save EMF (over limited spectral range).
% -------------------------------------------------------------------------
if params.save.emf == 1
    file_name_emf = [params.name,'.emf'];
    temp = figure();
    plot(freq1(display_indices)/1.0e9,10*log10(xf_display(display_indices)/1.0e-3),'LineWidth',2,'Color',[0 0 0]);
    axis off;
    print(temp,'-dmeta','-r300',file_name_emf);
    close(temp);
end
% Save the displayed spectrum in an EMF file (without axes)


% -------------------------------------------------------------------------
% Save ASCII (over limited spectral range).
% -------------------------------------------------------------------------
savearray = [(freq1(display_indices))'/1.09,(10*log10(xf_display(display_indices)/1.0e-3))'];
% Create array with the quantities to be saved.

if params.save.txt == 1
    file_name_ascii = [params.name,'.dat'];    
    fid = fopen(file_name_ascii,'w');
    fprintf(fid,'%s\t%s\n','frequency','power');
    fprintf(fid,'%s\t%s\n','(GHz)','(dBm)');
    fprintf(fid,'%g\t%g\n',savearray');
    fclose(fid);
end  
% Save the data of the displayed spectrum in a text file

end