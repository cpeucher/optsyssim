function char_opt_constellation(sig,params)
% Constellation diagram for optical signals
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function displays the constellation diagram of an optical signal
% (i.e. Im(E) vs Re(E) where E is the electric field sampled once per
% symbol).
% Rather use the dsp_constellation functions for digital constellation.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_constellation.type = 'scatter';%'transitions';%'scatter_transitions';
% params_constellation.pol = 'x';% 'y','both';
% params_constellation.nsamples_per_symbol = nsamples_per_symbol;
% params_constellation.sampling_time = nsamples_per_symbol/2;
% params_constellation.normalisation = 'mean';%'max';
% params_constellation.save.emf = 0;
% % params_constellation.line_color = 'b';
% % params_constellation.line_width = 3;
% % params_constellation.marker_type = 'o';
% % params_constellation.marker_color = 'r';
% % params_constellation.marker_size = 36;
% % params_constellation.plot_circular_grid = 1;
% % params_constellation.plot_axes = 0;
% params_constellation.name = 'Constellation';
% char_opt_constellation(sig,params_constellation); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               optical signal to characterise 
%                       [optical signal structure]
%
% params            optical constellation visualisation parameters
%                       [structure]
%
%                       params.type
%                           type of constellation to represent [string]
%
%                           params.type = 'scatter'
%                               only the samples at a specific sampling 
%                               time are represented
%                               One sample per symbol.
%
%                           params.type = 'transitions'
%                               transitions from symbol to symbol are 
%                               represented
%
%                           params.type = 'scatter_transitions'
%                               both samples and transitions are 
%                               represented.
%
%                           params.pol
%                               polarisation of the signal that is 
%                               represented [string]
%
%                               params.pol = 'x';
%                               params.pol = 'y';
%                               params.pol = 'both'; 
%
%                           params.nsamples_per_symbol
%                               number of samples per symbol in the signal 
%                               sig [integer scalar]
%
%                           params.sampling_time
%                               index of the sample that is represented in
%                               the 'scatter' mode [integer scalar]
%                               This index is relative to the symbol slot.
%
%                           params.normalisation
%                               type of normalisation for scatter modes
%                               [string]
%
%                               params.normalisation = 'max': maximum power
%                               params.normalisation = 'mean': mean powe
%
%                               In transitions modes, the constellation is 
%                               normalised to the peak power of the signal.
%
%                           params.name
%                               title of the visualiser [string]
%
%                           params.save
%                               specifies whether the constellation
%                               diagram should be saved as EMF file [0/1]
%                               params.save.emf = 0,1;
%
%                           Optional graphic parameters 
%                           (default values are indicated within brakets)
%
%                           params.line_color
%                               color of the line in transition modes ['b']
%
%                           params.line_width
%                               width of the line in transition modes [3]
%
%                           params.marker_type
%                               type of marker in scatter modes ['o']
%
%                           params.marker_color
%                               color of the markers in scatter modes ['b']
%
%                           params.marker_size
%                               marker size in scatter modes [36]
%
%                           params.plot_circular_grid
%                               plot circular grid [1]
%
%                           params.plot_axes
%                               plot axes [0]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% None.              
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% dt                time samples separation, in s [real scalar]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global dt

if isfield(params,'line_color')
    graphic_param.line_color = params.line_color;
else
    graphic_param.line_color = 'b';
end
    
if isfield(params,'marker_size')
    graphic_param.marker_size = params.marker_size;
else
    graphic_param.marker_size = 36;
end

if isfield(params,'line_width')
    graphic_param.line_width = params.line_width;
else
    graphic_param.line_width = 3;
end

if isfield(params,'marker_color')
    graphic_param.marker_color = params.marker_color;
else
    graphic_param.marker_color = 'r';
end

if isfield(params,'marker_type')
    graphic_param.marker_type = params.marker_type;
else
    graphic_param.marker_type = 'o';
end

if isfield(params,'plot_circular_grid')
    graphic_param.plot_circular_grid = params.plot_circular_grid
else
    graphic_param.plot_circular_grid = 1;
end

if isfield(params,'plot_axes')
    graphic_param.plot_axes = params.plot_axes;
else
    graphic_param.plot_axes = 0;
end
% In case of optional graphic parameters, overwrite the default value(s)


if strcmp(params.pol,'x')
    field = sig.x;
elseif strcmp(params.pol,'y')
    field = sig.y;
elseif strcmp(params.pol,'both')
    field = sqrt(abs(sig.x).^2 + abs(sig.y).^2) * exp(-(angle(sig.x) + angle(sig.y))/2);
end
% Define field to plot depending on polarisation

pav = mean(abs(field).^2);
% Average power of the field to plot


if strcmp(params.type,'scatter')
    
    samples = field(params.sampling_time:params.nsamples_per_symbol:length(field));
    % Collect samples of the signal
    
    if strcmp(params.normalisation,'max')
        pnorm = max(samples.*conj(samples));
    elseif strcmp(params.normalisation,'mean')
        pnorm = mean(samples.*conj(samples));
    else
        disp('char_opt_constellation: normalisation type not defined in scatter mode.');
    end
    % Normalisation
    
    figure('Name',params.name)
    hold on;
    scatter(gca,real(samples)/sqrt(pnorm),imag(samples)/sqrt(pnorm),'Line',2,'Marker',graphic_param.marker_type,'MarkerEdgeColor',graphic_param.marker_color,'MarkerFaceColor',graphic_param.marker_color,'SizeData',graphic_param.marker_size);
    
elseif strcmp(params.type,'transitions')
    
    pnorm = max(abs(field).^2);
    % Normalisation is the peak power
    
    figure('Name',params.name);
    hold on;
    plot(real(field)/sqrt(pnorm),imag(field)/sqrt(pnorm),'LineWidth',graphic_param.line_width);
    % Plot transitions
    
elseif strcmp(params.type,'scatter_transitions')
    
    pnorm = char_opt_peak_power(sig,[0 (length(sig.x)-1)*dt]);
    % Normalisation is the peak power
    
    samples = field(params.sampling_time:params.nsamples_per_symbol:length(field));
    % Collect samples of the signal
    
    figure('Name',params.name);
    hold on;
    plot(real(field)/sqrt(pnorm),imag(field)/sqrt(pnorm),'LineWidth',graphic_param.line_width,'Color',graphic_param.line_color);
    % Plot transitions
    scatter(gca,real(samples)/sqrt(pnorm),imag(samples)/sqrt(pnorm),'Line',2,'Marker',graphic_param.marker_type,'MarkerEdgeColor',graphic_param.marker_color,'MarkerFaceColor',graphic_param.marker_color,'SizeData',graphic_param.marker_size);
    % Plot states
    
    
else
    disp('char_opt_constellation: visualisation type not implemented.');
end


if graphic_param.plot_circular_grid == 1
    char_opt_constellation_plot_grid();
end
% Plot grid

axis([-1.1 1.1 -1.1 1.1]);
% Set axes range
axis square;
% Make orthonormal basis

if graphic_param.plot_axes == 1
    xlabel('real(E)');
    ylabel('imag(E)');
else
    axis off;
end

set(gca,'ytick',[]);
set(gca,'xtick',[]);
% Remove tick lines

if params.save.emf == 1
    print(gcf,'-dmeta','-r300',[params.name,'.emf']);
end
% Export to emf file

end
% End of main function



function char_opt_constellation_plot_grid()
x = 0:2*pi/50:2*pi;
plot(sin(x),cos(x),'-k');
plot(0.75*sin(x),0.75*cos(x),':k');
plot(0.5*sin(x),0.5*cos(x),'-k');
plot(0.25*sin(x),0.25*cos(x),':k');
% Display circles on graph
y = [-1 1];
plot(y,zeros(length(y),1),':k');
plot(zeros(length(y),1),y,':k');
plot(sqrt(2)/2*y,sqrt(2)/2*y,':k');
plot(sqrt(2)/2*y,-sqrt(2)/2*y,':k');
% Display grid on graphs
end
% End of plot circular grid auxilliary function