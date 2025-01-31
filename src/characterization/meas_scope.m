function varargout = meas_scope(sig,params)
% Visualisation of optical and electrical signals
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This functions visualises electrical or optical signals in the time
% domain.
% It is to be used for quick inspection of signals, since there is little
% control on the plot properties.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_scope.visualisers = {'amplitude','power','phase','chirp'};
% params_scope.display_interval = [0 time_array(end)];
% params_scope.save.emf = 0;
% params_scope.name = 'Waveform';
% meas_scope(sig,params_scope); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               electrical or optical signal to be visualised 
%                       [complex vector]
%               
%                       For an optical signal structure, sig can be
%                           sig.x to represent the signal along -x
%                           sig.y to represent the signal along -y
% 
%                       In case a single signal is split along the 2
%                       polarisation axes, one should construct:
%                           sqrt(abs(sig.x).^+ abs(sig.y).^2).*exp(-1i*phi)
%                           phi = -(angle(sig.x) + angle(sig.y))/2
%
% params            scope parameters [parameters]
%
%                       params.visualisers
%                           visualisers to consider [celle array]
%
%                           params.visualisers is a subset of 
%                               {'amplitude','power','phase','chirp'}
%
%                           'amplitude' rather applies to electrical
%                               signals
% 
%                           'power', 'phase', 'chirp' apply to optical
%                           signals
% 
%                           The visualisers in the list can be in any
%                           order. They will be sorted according to
%                           {'amplitude','power','phase','chirp'}
% 
%                       params.display_interval
%                           time interval over which the signal will be 
%                           displayed [2-element vector] [start stop]
%
%                       params.name
%                           name of the figure [string]
%
%                           Will also be used as file name basis
%
%                       params.save
%                           saving information [structure]
%
%                               params.save.emf
%                                   switch determining whether to save the 
%                                   trace of the waveform in an .emf file,
%                                   without legend nor axes, in order to be
%                                   directly used in e.g. pptx 
%                                   presentations) [0/1] 
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% varargout         optional export of trace data [structure]
%
%                       varargout{1}.time_chirp
%                           time axis values for chirp display, in s
%                           [real vector]
%
%                       varargout{1}.chirp
%                           frequency chirp, in Hz [real vector]
%
%                       varargout{1}.time
%                           time axis for power and phase display, in s
%                           [real vector]
%
%                       varargout{1}.power
%                           power, in W [real vector]
%
%                       varargout{1}.phase
%                           phase, in radians [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% time_array        time samples, in s [real vector]
% 
% dt                    time samples separation, in s [real scalar]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global time_array
global dt

plot_options = {'amplitude','power','phase','chirp'};
% List of all possible plots

visualisers = plot_options(ismember(plot_options,params.visualisers));
% Sort the quantities ones wishes to plot according to the order of the
% list above

nplots = length(visualisers);
% Number of plots


params.display_interval = params.display_interval(:).';
% Force display_interval to a line vector

if params.display_interval(1) > params.display_interval(2)
    fliplr(params.display_interval);
end
% Ensure the display interval is in increasing value order

display_indices = time_array <= params.display_interval(2) & time_array >= params.display_interval(1);
% Indices of time_array corresponding to the display range

display_time = cell(1,2);
display_time{1} = time_array(display_indices);
% Corresponding time values

display_time_interval = [display_time{1}(1) display_time{1}(end)];
sig = sig(:).';
% Ensure sig is a line vector

nsamples = length(sig);
% Number of samples in the input signal.


display_waveform = cell(1,4);
% Create empty cell array in which the waveforms to be plotted will be
% saved

if any(strcmp(visualisers,'amplitude')) 
    display_waveform{1} = sig(display_indices);
end

if any(strcmp(visualisers,'power'))
    display_waveform{2} = abs(sig(display_indices)).^2;
end

if any(strcmp(visualisers,'phase')) || any(strcmp(visualisers,'chirp'))
    signal_phase = unwrap(-angle(sig));
    display_waveform{3}  = signal_phase(display_indices);
end
% Calculate the signal phase, if necessary

if any(strcmp(visualisers,'chirp'))
    signal_chirp = -num_diff_1d_pb(signal_phase,dt)/2/pi;
    % Extract the chirp from the phase
    signal_chirp = signal_chirp(2:nsamples-1);
    % Restrict the chirp values to ensure there is no discontinuity
    % problem due to the fact that the phase is not periodic
    time_array_chirp = time_array(2:nsamples-1);
    % Corresponding time range

    display_indices_chirp = time_array_chirp <= params.display_interval(2) & time_array_chirp >= params.display_interval(1);
    % Indices of time_array_chirp corresponding to the display range
    display_time{2} = time_array_chirp(display_indices_chirp);
    % Corresponding time values

    display_waveform{4} = signal_chirp(display_indices_chirp);
end
% Calculate the signal chirp, if necessary




% Plot visualisers:

hfig.fig = figure('Name',params.name);

for iplot = 1:nplots

    hfig.ax(iplot) = subplot(nplots,1,iplot);

    plot_visualiser(visualisers{iplot},display_time,display_waveform,display_time_interval)

end
% End of loop over plots


if params.save.emf

    file_name_emf = [params.name,'.emf'];

    hfig_emf = figure;

    for iplot = 1:nplots
        copyobj(hfig.ax(iplot), hfig_emf);
    end

    set(get(gcf,'Children'),'Visible','off');

    print(hfig_emf,'-dmeta','-r300',file_name_emf);

    close(hfig_emf)

end
% Save the traces as an emf file
% For this purpose we duplicate the figure, then turn all axes off, before
% saving to emf.


% Optional export of trace data:
varargout{1}.time_chirp = display_time{2};
varargout{1}.chirp = display_waveform{4};
varargout{1}.time = display_time{1};
varargout{1}.power = display_waveform{2};
varargout{1}.phase = display_waveform{3};

end
% End of main function




% -------------------------------------------------------------------------
% Plot trace, depending on visualiser type
% -------------------------------------------------------------------------
function plot_visualiser(type,display_time,display_waveform,display_time_interval)

switch type

    case 'amplitude'

        plot(display_time{1}/1.0e-12,display_waveform{1},'b');
        grid on;
        xlabel('time (ps)');
        ylabel('amplitude (a.u.)');
        xlim(display_time_interval/1.0e-12)
        ylim([1.1*min(display_waveform{1}) - 0.1*max(display_waveform{1}) - eps ...
            1.1*max(display_waveform{1}) - 0.1*min(display_waveform{1}) + eps])

    case 'power'

        plot(display_time{1}/1.0e-12,display_waveform{2}/1.0e-3,'b');
        grid on;
        xlabel('time (ps)');
        ylabel('power (mW)');
        xlim(display_time_interval/1.0e-12)
        ylim([0 max(display_waveform{2})*1.1]/1.0e-3)

    case 'phase'

        plot(display_time{1}/1.0e-12,display_waveform{3},'b');
        grid on;
        xlabel('time (ps)');
        ylabel('phase (rad)');
        xlim(display_time_interval/1.0e-12)
        ylim([1.1*min(display_waveform{3}) - 0.1*max(display_waveform{3}) - eps ...
            1.1*max(display_waveform{3}) - 0.1*min(display_waveform{3}) + eps])


    case 'chirp'

        plot(display_time{2}/1.0e-12,display_waveform{4}/1.0e9,'b');
        grid on;
        xlabel('time (ps)');
        ylabel('chirp (GHz)');
        xlim(display_time_interval/1.0e-12)
        ylim([1.1*min(display_waveform{4}) - 0.1*max(display_waveform{4}) - eps ...
            1.1*max(display_waveform{4}) - 0.1*min(display_waveform{4}) + eps]/1.0e9);


end
% End of switch over plot types

end
% End of plot function








%     ----------------------------------------------------------------------
%     Range limitation for display
%     ----------------------------------------------------------------------     
%     if max(display_waveform) == 0
%         display_power_ymax = 1;
%     else
%         display_power_ymax = 1.2*max(display_waveform);
%     end 
%     display_power_ymin = 0;
%     display_power_range = [display_power_ymin display_power_ymax]/1e-3;
%     Set visualisation range for the power display, in mW
% 
%     display_phase_range(1) = min(display_phase) - 0.1;  
%     display_phase_range(2) = max(display_phase) + 0.1;
%     Set visualisation range for the phase display, in rad
% 
%     display_time_range = params.display_interval/1.0e-12;
%     Set time visualisation range, in ps 
% 
% 
% 

% 
%     ----------------------------------------------------------------------
%     Save into text file
%     ----------------------------------------------------------------------
% 
%     if params.save.txt == 1
%         save_array = [time_array',signal_waveform',signal_phase'];
%         save_array_chirp = [time_chirp',signal_chirp'];    
%         Build arrays containing the information to be saved in text 
%         files       
% 
%         file_name_ascii_1 = [params.name,'.dat'];        
%         fid = fopen(file_name_ascii_1,'w');
%         fprintf(fid,'%s\t%s\t%s\n','time','power','phase');
%         fprintf(fid,'%s\t%s\t%s\n','(s)','(W)','(rad)');
%         fprintf(fid,'%g\t%g\t%g\n',save_array');
%         fclose(fid);
% 
%         file_name_ascii_2 = [params.name,'_chirp.dat'];
%         fid = fopen(file_name_ascii_2,'w');
%         fprintf(fid,'%s\t%s\n','time','chirp');
%         fprintf(fid,'%s\t%s\n','(s)','(Hz)');
%         fprintf(fid,'%g\t%g\n',save_array_chirp');
%         fclose(fid);      
%     end
%     Save into text files
% 
% 
% 
% 

% 
% 
%     ----------------------------------------------------------------------
%     Save into text file
%     ----------------------------------------------------------------------    
%     if params.save.txt == 1
%         save_array=[time_array',signal_waveform'];
%         Build array containing the information to be saved in the text
%         file
% 
%         file_name_ascii = [params.name,'.dat'];
%         fid = fopen(file_name_ascii,'w');
%         fprintf(fid,'%s\t%s\n','time','amplitude');
%         fprintf(fid,'%s\t%s\n','(s)','(a.u.)');
%         fprintf(fid,'%g\t%g\n',save_array');
%         fclose(fid);
%     end
%     Save into text file
% 
% 
% 
