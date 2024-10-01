function varargout = meas_scope(sig,params)
% Oscilloscope for visualisation of optical and electrical signals
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
% params_scope.type = 'opt';%'elec';
% params_scope.pol = 'x';%'y','both';
% params_scope.visualiser_type = 'power';%'power_phase';'power_chirp';'power_phase_chirp';
% params_scope.display_interval = [0 time_array(end)];
% params_scope.save.txt = 0;
% params_scope.save.emf = 0;
% params_scope.name = 'Waveform';
% meas_scope(sig,params_scope); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               electrical or optical signal to be visualised 
%                       [real vector or optical signal structure]
%
% params            scope parameters [parameters]
%
%                       params.type
%                           type of  signal [string]
%
%                           params.type = 'elec' for electrical signals
%                           params.type = 'opt' for optical signals
%
%                       params.pol
%                           polarisation of the signal that will be
%                           visualised, in case it is optical [string]
%
%                           params.pol = 'x';
%                           params.pol = 'y';
%                           params.pol = 'both';
%
%                       params.visualiser_type
%                           type of visualisation [string]
% 
%                           Applies to optical signals.
%
%                           params.visualiser_type = 'power'
%                           params.visualiser_type = 'power_phase'
%                           params.visualiser_type = 'power_chirp'
%                           params.visualiser_type = 'power_phase_chirp'
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
%                               params.save.txt
%                                   switch determining whether to save the 
%                                   waveform in a .dat text file [0/1] 
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


params.display_interval = params.display_interval(:).';
% Force display_interval to a line vector

if params.display_interval(1) > params.display_interval(2)
    fliplr(params.display_interval);
end
% Ensure the display interval is in increasing value order


%--------------------------------------------------------------------------
% Pre-processing of optical signals.
%--------------------------------------------------------------------------
if strcmp(params.type,'opt')
    % Process optical signals
    nsamples = length(sig.x);
    % Number of samples in the optical signal
    
    % First, prepare the quantities of interest depending on the
    % polarisation that the user wishes to plot   
    if strcmp(params.pol,'x')
        % One wants to visualise the contribution of the signal polarised
        % along -x
        signal_waveform = abs(sig.x).^2;
        % Intensity
        signal_phase = unwrap(-angle(sig.x));
        % Phase
    elseif strcmp(params.pol,'y')
        % One wants to visualise the contribution of the signal polarised
        % along -y
        signal_waveform = abs(sig.y).^2;
        % Intensity
        signal_phase = unwrap(-angle(sig.y));
        % Phase
    elseif strcmp(params.pol,'both')
        % The signal has an arbitrary state of polarisation.
        signal_waveform = abs(sig.x).^2 + abs(sig.y).^2;
        % Intensity
        signal_phase = unwrap(-(angle(sig.x) + angle(sig.y))/2);
        % Common phase
    end
    
    
    % Next, the phase is processed in order to calculate the chirp    
    signal_chirp = -num_diff_1d_pb(signal_phase,dt)/2/pi;
    % Extract the chirp from the phase
    signal_chirp = signal_chirp(2:nsamples-1);
    % Restrict the chirp values to ensure there is no discontinuity
    % problem due to the fact that the phase is not periodic
    time_chirp = time_array(2:nsamples-1);
    % Corresponding time range
    
    
    % Next, the range of the various quantities is restricted to match the
    % desired visualisation range
    range_indices = find(time_array >= params.display_interval(1) & time_array <= params.display_interval(2));
    % Find the indices of the time array corresponding to the specified time
    % interval
    range_indices_chirp = find(time_chirp >= params.display_interval(1) & time_chirp <= params.display_interval(2));
    % Same thing for the time array corresponding to the chirp
    
    display_time = time_array(range_indices);
    display_waveform = signal_waveform(range_indices);
    display_phase = signal_phase(range_indices);
    % Restrict the intensity and phase to the visualisation range
    
    display_time_chirp = time_chirp(range_indices_chirp);
    display_chirp = signal_chirp(range_indices_chirp);
    % Restrict the chirp to the visualisation range
    
    varargout{1}.time_chirp = display_time_chirp;
    varargout{1}.chirp = display_chirp;
    varargout{1}.time = display_time;
    varargout{1}.power = display_waveform;
    varargout{1}.phase = display_phase;
    % Save traces in optional variable argument structure
    
    %----------------------------------------------------------------------
    % Range limitation for display
    %----------------------------------------------------------------------     
    if max(display_waveform) == 0
        display_power_ymax = 1;
    else
        display_power_ymax = 1.2*max(display_waveform);
    end 
    display_power_ymin = 0;
    display_power_range = [display_power_ymin display_power_ymax]/1e-3;
    % Set visualisation range for the power display, in mW
    
    display_phase_range(1) = min(display_phase) - 0.1;  
    display_phase_range(2) = max(display_phase) + 0.1;
    % Set visualisation range for the phase display, in rad
    
    display_time_range = params.display_interval/1.0e-12;
    % Set time visualisation range, in ps 
    
    
    
    %----------------------------------------------------------------------
    % Power visualisation
    %----------------------------------------------------------------------       
    if strcmp(params.visualiser_type,'power')
        % 'power' display
        
        figure('Name',params.name)
        plot(display_time/1e-12,display_waveform/1.0e-3,'Color',[0 0 1]);
        xlim(display_time_range);
        ylim(display_power_range);
        grid on;
        xlabel('time (ps)');
        ylabel('power (mW)');
        % Standard screen display
        
        if params.save.emf == 1
            file_name_emf = [params.name,'.emf'];
            temp = figure();
            plot(display_time/1e-12,display_waveform/1.0e-3,'LineWidth',2,'Color',[0 0 0]);
            xlim(display_time_range);
            ylim(display_power_range);
            axis off;
            print(temp,'-dmeta','-r300',file_name_emf);
            close(temp);
        end
        % Save traces only into emf file
        
        
        
    %----------------------------------------------------------------------
    % Power-phase visualisation
    %----------------------------------------------------------------------        
    elseif strcmp(params.visualiser_type,'power_phase')
        % 'power_phase' display
        
        figure('Name',params.name)
        subplot(2,1,1)
        plot(display_time/1e-12,display_waveform/1.0e-3,'Color',[0 0 1]);
        xlim(display_time_range);
        ylim(display_power_range);
        grid on;
        xlabel('time (ps)');
        ylabel('power (mW)');
        subplot(2,1,2)
        plot(display_time/1e-12,display_phase,'Color',[0 0 1]);        
        xlim(display_time_range);
        ylim(display_phase_range);     
        grid on;
        xlabel('time (ps)');
        ylabel('phase (rad)');
        % Standard screen display
        
        
        if params.save.emf == 1
            file_name_emf = [params.name,'.emf'];
            temp = figure();
            subplot(2,1,1)
            plot(display_time/1e-12,display_waveform/1.0e-3,'Color',[0 0 1]);
            xlim(display_time_range);
            ylim(display_power_range);
            axis off;
            subplot(2,1,2)
            plot(display_time/1e-12,display_phase,'Color',[0 0 1]);
            xlim(display_time_range);
            ylim(display_phase_range);  
            axis off;
            print(temp,'-dmeta','-r300',file_name_emf);
            close(temp);
        end
        % Save traces only into emf file
        
    %----------------------------------------------------------------------
    % Power-chirp visualisation
    %----------------------------------------------------------------------         
    elseif strcmp(params.visualiser_type,'power_chirp')
        % 'power_chirp' display
        
        figure('Name',params.name)
        subplot(2,1,1)
        plot(display_time/1e-12,display_waveform/1.0e-3,'Color',[0 0 1]);
        xlim(display_time_range);
        ylim(display_power_range);
        grid on;
        xlabel('time (ps)');
        ylabel('power (mW)');
        subplot(2,1,2)
        plot(display_time_chirp/1e-12,display_chirp/1.0e9,'Color',[0 0 1]);
        xlim(display_time_range);
        grid on;
        xlabel('time (ps)');
        ylabel('chirp (GHz)');
        % Standard screen display
        
        if params.save.emf == 1
            file_name_emf=[params.name,'.emf'];
            temp = figure();
            subplot(2,1,1)
            plot(display_time/1e-12,display_waveform/1.0e-3,'Color',[0 0 0]);
            xlim(display_time_range);
            ylim(display_power_range);
            axis off;
            subplot(2,1,2)
            plot(display_time_chirp/1e-12,display_chirp/1.0e9,'Color',[0 0 0]);
            xlim(display_time_range);
            axis off;
            print(temp,'-dmeta','-r300',file_name_emf);
            close(temp);
        end
        % Save traces only into emf file
        
    %----------------------------------------------------------------------
    % power-phase-chirp visualisation
    %----------------------------------------------------------------------        
    elseif strcmp(params.visualiser_type,'power_phase_chirp')
        % 'power_phase_chirp' display
        
        figure('Name',params.name)
        subplot(3,1,1)
        plot(display_time/1e-12,display_waveform/1.0e-3,'Color',[0 0 1]);
        xlim(display_time_range);
        ylim(display_power_range);
        grid on;
        xlabel('time (ps)');
        ylabel('power (mW)');
        subplot(3,1,2)
        plot(display_time/1e-12,display_phase,'Color',[0 0 1]);
        xlim(display_time_range);
        ylim(display_phase_range);
        grid on;
        xlabel('time (ps)');
        ylabel('phase (rad)');
        subplot(3,1,3)
        plot(display_time_chirp/1e-12,display_chirp/1.0e9,'Color',[0 0 1]);
        xlim(display_time_range);
        grid on;
        xlabel('time (ps)');
        ylabel('chirp (GHz)');
        % Standard screen display
        
        
        if params.save.emf == 1
            file_name_emf = [params.name,'.emf'];
            temp = figure();
            subplot(3,1,1)
            plot(display_time/1e-12,display_waveform/1.0e-3,'Color',[0 0 1]);
            xlim(display_time_range);
            ylim(display_power_range);
            axis off;
            subplot(3,1,2)
            plot(display_time/1e-12,display_phase,'Color',[0 0 1]);
            xlim(display_time_range);
            ylim(display_phase_range);
            axis off;
            subplot(3,1,3)
            plot(display_time_chirp/1e-12,display_chirp/1.0e9,'Color',[0 0 1]);
            xlim(display_time_range);
            axis off;
            print(temp,'-dmeta','-r300',file_name_emf);
            close(temp);
        end   
        % Save traces only into emf file
        
        
    else
        error('\n\n%s\n','meas_scope: visualisation not implemented for optical signals.')
    end
    % End of choice of visualisation type
    
    
    %----------------------------------------------------------------------
    % Save into text file
    %----------------------------------------------------------------------
    
    if params.save.txt == 1
        save_array = [time_array',signal_waveform',signal_phase'];
        save_array_chirp = [time_chirp',signal_chirp'];    
        % Build arrays containing the information to be saved in text 
        % files       
        
        file_name_ascii_1 = [params.name,'.dat'];        
        fid = fopen(file_name_ascii_1,'w');
        fprintf(fid,'%s\t%s\t%s\n','time','power','phase');
        fprintf(fid,'%s\t%s\t%s\n','(s)','(W)','(rad)');
        fprintf(fid,'%g\t%g\t%g\n',save_array');
        fclose(fid);
        
        file_name_ascii_2 = [params.name,'_chirp.dat'];
        fid = fopen(file_name_ascii_2,'w');
        fprintf(fid,'%s\t%s\n','time','chirp');
        fprintf(fid,'%s\t%s\n','(s)','(Hz)');
        fprintf(fid,'%g\t%g\n',save_array_chirp');
        fclose(fid);      
    end
    % Save into text files
    
    
    
    
%--------------------------------------------------------------------------
% Electrical signals
%--------------------------------------------------------------------------    
elseif strcmp(params.type,'elec')
    % Process electrical signals
    
    %----------------------------------------------------------------------
    % Plot electrical waveform
    %----------------------------------------------------------------------
    signal_waveform = sig;
    % If the signal is electrical then the voltage/current is plotted
    
    range_indices = find(time_array >= params.display_interval(1) & time_array <= params.display_interval(2));
    % Find the indices of the time array corresponding to the specified time
    % interval
    
    display_time = time_array(range_indices);
    display_waveform = signal_waveform(range_indices);
    % Extract the samples corresponding to the desired time interval
    
    figure('Name',params.name)
    plot(display_time/1e-12,display_waveform,'Color',[0 0 1]);
    grid on;
    xlabel('time (ps)');
    ylabel('amplitude (a.u.)');
    % Standard screen display    

    
    if params.save.emf == 1
        file_name_emf = [params.name,'.emf'];
        temp = figure();
        plot(display_time/1e-12,display_waveform,'LineWidth',2,'Color',[0 0 0]);
        axis off;
        print(temp,'-dmeta','-r300',file_name_emf);
        close(temp);
    end
    % Save traces only into emf file
    
    


    %----------------------------------------------------------------------
    % Save into text file
    %----------------------------------------------------------------------    
    if params.save.txt == 1
        save_array=[time_array',signal_waveform'];
        % Build array containing the information to be saved in the text
        % file
        
        file_name_ascii = [params.name,'.dat'];
        fid = fopen(file_name_ascii,'w');
        fprintf(fid,'%s\t%s\n','time','amplitude');
        fprintf(fid,'%s\t%s\n','(s)','(a.u.)');
        fprintf(fid,'%g\t%g\n',save_array');
        fclose(fid);
    end
    % Save into text file
    
    
    
end
% End of choice of optical or electrical signal

end