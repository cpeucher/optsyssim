function varargout = meas_eye(sig,params)
% Visualisation of eye diagram of optical or electrical signals
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function produces an eye diagram 
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_eye.pol = 'x';%'y','both';
% params_eye.neyes = 2;
% params_eye.nsamples_per_symbol = nsamples_per_symbol;
% params_eye.save.txt = 0;
% params_eye.save.emf = 0;
% params_eye.save.jpg = 0;
% params_eye.save.vertical_scale_type = 'auto';%'fixed';
% params_eye.save.vertical_scale_max = 1.0e-3;
% params_eye.save.display_baseline = 1;
% params_eye.colour_grade = 0;
% params_eye.name = 'Eye diagram';
% meas_eye(sig,params_eye);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               signal to be visualised
%                       [real vector]
%
% params            eye diagram parameters [structure]
%
%                        params.pol
%                           polarisation of the signal that will 
%                           be visualised, in case it is optical [string]
%
%                           params.pol = 'x'
%                           params.pol = 'y' 
%                           params.pol = 'both'
%
%                       params.neyes
%                           number of eyes in the eye diagram 
%                           [integer scalar]
%
%                       params.nsamples_per_symbol
%                           number of samples per symbol of the the signal
%                           to be visualised [integer scalar]
%
%                       params.name
%                           name of the figure [string]
%                           Will also be used as file name basis.
%
%                       params.save
%                           saving information [structure]
%
%                               params.save.txt
%                                   switch to save the eye diagram in 
%                                   .dat text file consisting of the time
%                                   axis + one column per trace in the eye
%                                   diagram [0/1]
%
%                               params.save.emf
%                                   switch to save the eye diagram in 
%                                   .emf file (without legend nor axes, in 
%                                   order to be directly used in e.g. ppt 
%                                   presentation) [0/1]   
%
%                               params.save.vertical_scale_type
%                                   specify whether the vertical scale is 
%                                   adapted to the signal value, or whether 
%                                   it is kept at a fixed value, allowing 
%                                   for comparison between different eye
%                                   diagrams [string]
%
%                                   params.save.vertical_scale_type = 'auto';
%                                   params.save.vertical_scale_type = 'fixed';
%
%                               params.save.vertical_scale_max
%                                   maximum value of the eye diagram in 
%                                   case 
%                                   params.save.vertical_scale_type='fixed'
%                                   [real scalar]
%
%                               params.save.display_baseline
%                                   specifies whether to display the 
%                                   baseline in the saved eye diagram [0/1]
%
%                                   params.save.display_baseline = 1;
%                                   params.save.display_baseline = 0;
%
%                               params.colour_grade
%                                   specify whether a colour grade eye
%                                   diagram will be visualised and saved
%                                   [0/1]
%
%                                   params.colour_grade = 1;
%                                   params.colour_grade = 0;  
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% varargout(1)      time base, in s [real vector]
%
% varargout(2)      traces [real matrix]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% dt                 time samples separation, in s [real scalar]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global dt


if isoptical(sig)
    % Test if the signal is an optical signal

    if strcmp(params.pol,'x')
        waveform = abs(sig.x).^2;

    elseif strcmp(params.pol,'y')
        waveform = abs(sig.y).^2;

    elseif strcmp(params.pol,'both')
        waveform = abs(sig.x).^2 + abs(sig.y).^2;

    end        
    % If the signal is optical, the quantity to plot is the power

else
    % Else it is an electrical signal

    waveform = sig;
    % If the signal is electrical then the voltage/current is plotted.    
end


nsymbols = length(waveform)/params.nsamples_per_symbol;
% Number of symbols
% If global parameters are set properly, this should be an integer.

nsamples_per_trace = params.nsamples_per_symbol*params.neyes;
% Number of samples per trace

eye_traces = reshape(waveform,nsamples_per_trace,nsymbols/params.neyes);
% Traces of the eye diagram
eye_time = (0:nsamples_per_trace - 1)*dt;
% Time base

varargout(1) = {eye_time};
varargout(2) = {eye_traces};
% Save to optional output parameters


% -------------------------------------------------------------------------
% Display
% -------------------------------------------------------------------------
min_max = abs(max(waveform) - min(waveform));
if min_max < eps
    display_min = 0;
    display_max = max(waveform)*1.1;
else
    display_delta = abs(max(waveform) - min(waveform))*0.1;
    display_min = min(waveform) - display_delta;
    display_max = max(waveform) + display_delta;
end
% Determine vertical range of the display


hfig_line_eye = figure('Name',params.name);
plot(eye_time./1e-12,eye_traces,'Color',[0 0 1],'LineWidth',2);
axis([0 params.neyes*params.nsamples_per_symbol*dt/1.0e-12 display_min display_max])
xlabel('time (ps)');
ylabel('amplitude / power (a.u.)');

if params.save.jpg == 1
    file_name_jpg = [params.name,'.jpg'];    
    print(hfig_line_eye,'-djpeg','-r300',file_name_jpg); 
end
    
if params.save.txt == 1
    eye_diagram_txt = [eye_time',eye_traces];
    file_name_txt = [params.name,'.dat'];
    save(file_name_txt,'eye_diagram_txt','-ascii','-double','-tabs')    
end

if params.save.emf == 1
    file_name_emf = [params.name,'.emf'];
    temp = figure();
    plot(eye_time./1e-12,eye_traces,'LineWidth',2,'Color',[0 0 0]);
    if min(waveform) < 0
        vertical_axis_min = min(waveform);
    else
        vertical_axis_min = -1.0e-5;
    end
    if strcmp(params.save.vertical_scale_type,'auto')
        axis([eye_time(1)/1.0e-12 eye_time(length(eye_time))/1.0e-12 vertical_axis_min max(waveform)]);
    elseif strcmp(params.save.vertical_scale_type,'fixed')
        axis([eye_time(1)/1.0e-12 eye_time(length(eye_time))/1.0e-12 vertical_axis_min params.save.vertical_scale_max]);
    else
        error('meas_eye: vertical scale type not defined.');
    end
    if params.save.display_baseline == 1
        hold on;
        plot(eye_time./1e-12,zeros(1,length(eye_time)),'LineWidth',1,'Color',[0 0 0],'LineStyle','--');
    end
    axis off;
    print(temp,'-dmeta','-r300',file_name_emf);
    close(temp);
end


% -------------------------------------------------------------------------
% Display and save colour grade eye diagram
% -------------------------------------------------------------------------

% TODO Clean up this mess and use histcounts2 


if params.colour_grade == 1
    % Check if a colour grade eye diagram is expected.
    
    amp_vert = max(waveform) - min(waveform);
    vertical_range_max = max(waveform) + 0.1*amp_vert;
    vertical_range_min = min(waveform) - 0.1*amp_vert;
    % Determine minimum and maximum values for the vertical scale.

    n_horizontal_pixels = nsamples_per_trace - 1;
    n_vertical_pixels = 100;
    % Number of pixels; by default the number of horizontal coordinates is
    % equal to the number of samples per trace.

    border_vert = vertical_range_min + (0:n_vertical_pixels)*(vertical_range_max - vertical_range_min)/n_vertical_pixels;
    % Array contains pixels vertical boundaries; NPixelsVert+1 elements.

    border_hori = 1 + (0:n_horizontal_pixels)*(nsamples_per_trace - 1)/n_horizontal_pixels;
    % Array containing pixels horizontal boundaries, in terms of the 
    % indices i of eyetime(i) or eyetrace(i,j); NPixelsHor+1 elements.
    % Observe that the elements of BorderHor are not necessarily integer, 
    % unless the number of horizontal pixels is equal to the number of 
    % samples per trace minus 1.

    % to remember:
    % eyetraces(:,i) contains all the samples of trace number i (of length
    % nsamples_per_trace
    % eyetraces(j,:) contains the collection of samples # j of all traces
    % (length is Nbits/params.neyes)

    count_matrix = zeros(length(border_vert),n_horizontal_pixels);
    % Initialise the matrix in which the values of the pixels will be 
    % stored.
    
    sample_range = (1:nsamples_per_trace);
    % Array containing the full range of horizontal samples in the eye 
    % diagram.
    for i=1:n_horizontal_pixels
        % For each horizontal bin.
        if i < n_horizontal_pixels
            sample_range_bin = find(sample_range >= border_hori(i) & sample_range < border_hori(i+1));
        elseif i == n_horizontal_pixels
            sample_range_bin = find(sample_range >= border_hori(i) & sample_range <= border_hori(i+1));            
        end
        
        for k = 1:length(sample_range_bin)
            % For each time sample in the bin.
            count_matrix(:,i) = count_matrix(:,i) + histc(eye_traces(sample_range_bin(k),:),border_vert)';
            % Count the number of values in each vertical bin for that
            % particular time sample.
        end
    end

    hfig_colour_grade_eye = figure('Name',[params.name ' (Colour grade)']);
    image([eye_time(1)/1.0e-12,nsamples_per_trace*dt/1.0e-12],[vertical_range_min, vertical_range_max],count_matrix);
    colormap(hot);
    set(gca,'YDir','normal');
    set(get(gca,'XLabel'),'String', 'time (ps)');
    set(get(gca,'YLabel'),'String', 'amplitude / power (a.u.)');
    % Visualise colour grade mode eye diagram
    

    
    if params.save.jpg == 1
        file_name_jpg = [params.Name,'_colour_grade','.jpg'];    
        print(hfig_colour_grade_eye,'-djpeg','-r300',file_name_jpg);
    end
    
    if params.save.emf == 1    
        %set(ColourGradeEye,'Units','centimeters','Position',[10 8 11 9]);
        set(hfig_colour_grade_eye,'PaperPositionMode','auto');
        file_name_emf=[params.Name,'_colour_grade','.emf'];    
        print(hfig_colour_grade_eye,'-dmeta','-r300',file_name_emf);
    end
        
end
% End of colour grade eye diagram plot.


end
