% -------------------------------------------------------------------------
% Test of opt_tf_fbg function
%
% Calculates different FBG transfer functions, including uniform, apodized
% and chirped.
%
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Preparation
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% ------------------------------------------------------------------------- 
% Clean up
% ------------------------------------------------------------------------- 
clear all
close all

% ------------------------------------------------------------------------- 
% Specify display format
% ------------------------------------------------------------------------- 
format long

% ------------------------------------------------------------------------- 
% Dock figures
% ------------------------------------------------------------------------- 
% set(0,'DefaultFigureWindowStyle','docked');
set(0,'DefaultFigureWindowStyle','normal');

% ------------------------------------------------------------------------- 
% Reinitialise the random number generator for reproducibility of the
% results
% ------------------------------------------------------------------------- 
stream = RandStream.getGlobalStream;
reset(stream);

% -------------------------------------------------------------------------
% Strings for filenames
% -------------------------------------------------------------------------
file_name_core_figure = strrep(mfilename,'run','fig');
file_name_core_data = strrep(mfilename,'run','data');
time_stamp = string(datetime('now','TimeZone','Z','Format','yyyyMMdd''T''HHmmss''Z'));

%--------------------------------------------------------------------------
% Switches
%--------------------------------------------------------------------------
do_debug = 1;
do_print = 0;
do_add_figsize_to_filename = 1;
margin_figure = 0;

% -------------------------------------------------------------------------
% Figure settings
% -------------------------------------------------------------------------
fig.font_name = 'times';
fig.interpreter = 'latex';
fig.font_size = 18;


mred = [178 0 24]/255;

% line_color = linspecer(6,'qualitative');
% Use plot(x,y,'Color',line_color(1,:));


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Global simulation parameters
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% ------------------------------------------------------------------------- 
% Define global variables
% ------------------------------------------------------------------------- 
global reference_frequency 
global frequency_array 
global time_array
global dt
global df 
global CONSTANT
% global space_grid

% -------------------------------------------------------------------------
% Set global simulation parameters
% -------------------------------------------------------------------------
reference_frequency = 193.1e12;
nsamples_per_symbol = 128;
nsymbols = 128;
symbol_rate = 40e9;

nsamples = nsamples_per_symbol*nsymbols;
sample_rate = nsamples_per_symbol*symbol_rate;
[time_array,dt,frequency_array,df] = core_create_time_axis(nsamples_per_symbol,nsymbols,symbol_rate);        
        
  
% reference_frequency = 193.1e12;
% df = 10e6;
% nsamples = 2^14;
% 
% sample_rate = nsamples*df;
% dt = 1/sample_rate;
% time_array = (0:nsamples-1)*dt;
% frequency_array = (-nsamples/2:nsamples/2-1)*df;

% -------------------------------------------------------------------------
% Space grid
% -------------------------------------------------------------------------
% xrange = [-15e-6, 15e-6];  
% yrange = [-15e-6, 15e-6]; 
% nxpoints = 2001;
% nypoints = 2001;
% space_grid = create_space_grid(xrange,yrange,nxpoints,nypoints);
% [space_grid.THETA,space_grid.RHO] = cart2pol(space_grid.X,space_grid.Y); 

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Startup routines
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% ------------------------------------------------------------------------- 
% Load physical constants
% ------------------------------------------------------------------------- 
CONSTANT = core_load_constants();
% Load essential physical constants

% ------------------------------------------------------------------------- 
% Start time 
% ------------------------------------------------------------------------- 
start_time = datetime("now");
fprintf('\n\n%s%s\n\n','Simulation started on ',start_time);


% ------------------------------------------------------------------------- 
% ------------------------------------------------------------------------- 
% ------------------------------------------------------------------------- 
% ------------------------------------------------------------------------- 
% Now we are ready to implement the system
% -------------------------------------------------------------------------        
% ------------------------------------------------------------------------- 
% ------------------------------------------------------------------------- 
% ------------------------------------------------------------------------- 



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------
% ------------------------------------------------------------------------- 
% -------------------------------------------------------------------------

%--------------------------------------------------------------------------
% WARNING: global parameters will be reset. Not to use for "real" system
% simulations.
% The goal here is just to visualise the fibre gratings transfer functions
% and to give the user full control of the wavelength range.
% This is motivated by the implementation of the calc_tf_dispersion
% function, which is currently designed to be used within the full
% simulation templates.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Wavelength range
%--------------------------------------------------------------------------
lambda_start = 1548e-9;
lambda_stop = 1552e-9;
lambda_step = 0.001e-9;

%--------------------------------------------------------------------------
% Grating definitions
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Uniform grating \kappa L= 2
%--------------------------------------------------------------------------
params_fbg1.length = 0.01;          % Grating length, in m.
params_fbg1.centre_frequency = CONSTANT.c/1550e-9;  % Centre frequency of the grating, in Hz.
params_fbg1.n0 = 1.45;              % Effective index of the fibre.
params_fbg1.dn = 1.0e-4;            % Refractive index change.
params_fbg1.m = 1;                  % Control of average index change.
params_fbg1.chirp_rate = 0;         % Linear chirp rate, in 1/m.
params_fbg1.apodisation.type = 'uniform';     
params_fbg1.apodisation.eta = 0;
params_fbg1.apodisation.fwhm = 0;
params_fbg1.apodisation.profile = 0;
numparams_fbg1.nsections = 1;       % Number of uniform grating sections.
params_fbg1.phase_shift = zeros(numparams_fbg1.nsections);% Position of phase shifts.


%--------------------------------------------------------------------------
% Uniform grating \kappa L= 8
%--------------------------------------------------------------------------
params_fbg2.length = 0.01;          % Grating length, in m.
params_fbg2.centre_frequency = CONSTANT.c/1550e-9;  % Centre frequency of the grating, in Hz.
params_fbg2.n0 = 1.45;             % Effective index of the fibre.
params_fbg2.dn = 4.0e-4;            % Refractive index change.
params_fbg2.m = 1;                  % Control of average index change.
params_fbg2.chirp_rate = 0;          % Linear chirp rate, in 1/m.
params_fbg2.apodisation.type = 'uniform';     
params_fbg2.apodisation.eta = 0;
params_fbg2.apodisation.fwhm = 0;
params_fbg2.apodisation.profile = 0;
numparams_fbg2.nsections = 1000;        % Number of uniform grating sections.
params_fbg2.phase_shift = zeros(numparams_fbg2.nsections);% Position of phase shifts.


%--------------------------------------------------------------------------
% Gaussian grating \kappa FWHM= 2
%--------------------------------------------------------------------------
params_fbg3.length = 0.03;          % Grating length, in m.
params_fbg3.centre_frequency = CONSTANT.c/1550e-9;  % Centre frequency of the grating, in Hz.
params_fbg3.n0 = 1.45;              % Effective index of the fibre.
params_fbg3.dn = 1.0e-4;            % Refractive index change.
params_fbg3.m = 1;                  % Control of average index change.
params_fbg3.chirp_rate = 0;          % Linear chirp rate, in 1/m.
params_fbg3.apodisation.type = 'gaussian';     
params_fbg3.apodisation.eta = 0;
params_fbg3.apodisation.fwhm = 0.01;
params_fbg3.apodisation.profile = 0;
numparams_fbg3.nsections = 1000;        % Number of uniform grating sections.
params_fbg3.phase_shift = zeros(numparams_fbg3.nsections);% Position of phase shifts.

%--------------------------------------------------------------------------
% Gaussian grating \kappa FWHM= 8
%--------------------------------------------------------------------------
params_fbg4.length = 0.03;          % Grating length, in m.
params_fbg4.centre_frequency = CONSTANT.c/1550e-9;  % Centre frequency of the grating, in Hz.
params_fbg4.n0 = 1.45;              % Effective index of the fibre.
params_fbg4.dn = 4.0e-4;            % Refractive index change.
params_fbg4.m = 1;                  % Control of average index change.
params_fbg4.chirp_rate = 0;          % Linear chirp rate, in 1/m.
params_fbg4.apodisation.type = 'gaussian';     
params_fbg4.apodisation.eta = 0;
params_fbg4.apodisation.fwhm = 0.01;
params_fbg4.apodisation.profile = 0;
numparams_fbg4.nsections = 1000;        % Number of uniform grating sections.
params_fbg4.phase_shift = zeros(numparams_fbg4.nsections);% Position of phase shifts.

%--------------------------------------------------------------------------
% Uniform chirped grating \kappa L= 20
%--------------------------------------------------------------------------
params_fbg5.length = 0.10;          % Grating length, in m.
params_fbg5.centre_frequency = CONSTANT.c/1550e-9;  % Centre frequency of the grating, in Hz.
params_fbg5.n0 = 1.45;              % Effective index of the fibre.
params_fbg5.dn = 1.0e-4;            % Refractive index change.
params_fbg5.m = 0;                  % Control of average index change.
params_fbg5.chirp_rate = 6.9e-9;          % Linear chirp rate, in 1/m.
params_fbg5.apodisation.type = 'uniform';     
params_fbg5.apodisation.eta = 0;
params_fbg5.apodisation.fwhm = 0;
params_fbg5.apodisation.profile = 0;
numparams_fbg5.nsections = 1000;        % Number of uniform grating sections.
params_fbg5.phase_shift = zeros(numparams_fbg5.nsections);% Position of phase shifts.


%--------------------------------------------------------------------------
% Quarter cosine chirped grating \kappa L= 20
%--------------------------------------------------------------------------
params_fbg6.length = 0.10;          % Grating length, in m.
params_fbg6.centre_frequency = CONSTANT.c/1550e-9;  % Centre frequency of the grating, in Hz.
params_fbg6.n0 = 1.45;              % Effective index of the fibre.
params_fbg6.dn = 1.0e-4;            % Refractive index change.
params_fbg6.m = 0;                  % Control of average index change.
params_fbg6.chirp_rate = -6.9e-9;          % Linear chirp rate, in 1/m.
params_fbg6.apodisation.type = 'quarter_cosine';     
params_fbg6.apodisation.eta = 0.1;
params_fbg6.apodisation.fwhm = 0;
params_fbg6.apodisation.profile = 0;
numparams_fbg6.nsections = 1000;        % Number of uniform grating sections.
params_fbg6.phase_shift = zeros(numparams_fbg6.nsections);% Position of phase shifts.



        
        
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% HERE WE GO
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------   

% First we redefine some of the global parameters.
fmin = CONSTANT.c/lambda_stop;
% Minimum frequency.
fmax = CONSTANT.c/lambda_start;
% Maximum frequency.
lambda_mean = (lambda_start+lambda_stop)/2;
% Centre wavelength

df = CONSTANT.c/lambda_mean.^2*lambda_step;
% Frequency step.
reference_frequency = CONSTANT.c/lambda_mean;
% Centre frequency.
freq = [fmin:df:fmax];
% Absolute frequency array.
frequency_array = freq - reference_frequency;
% Relative frequency array.


%--------------------------------------------------------------------------
% FBG1
%--------------------------------------------------------------------------
tf_fbg1 = opt_tf_fbg(freq,params_fbg1,numparams_fbg1);
% Calculate grating transfer function.

save_dispersion.status = 0;
save_dispersion.file_name = 'filter_tf.dat';
[phase_fbg1,freq_group_delay_fbg1,group_delay_fbg1,freq_dispersion_fbg1,dispersion_fbg1] = extract_dispersion_from_tf(frequency_array,tf_fbg1.transmission,'si',save_dispersion);
% Extract dispersion from transfer function.
% Calculate group delay and dispersion.
% We use the transfer function in transmission to avoid phase
% discontinuities. OK for symmetric apodisation profiles.


%--------------------------------------------------------------------------
% FBG2
%--------------------------------------------------------------------------
tf_fbg2 = opt_tf_fbg(freq,params_fbg2,numparams_fbg2);
% Calculate grating transfer function.


save_dispersion.status = 0;
save_dispersion.file_name = 'filter_tf.dat';
[phase_fbg2,freq_group_delay_fbg2,group_delay_fbg2,freq_dispersion_fbg2,dispersion_fbg2] = extract_dispersion_from_tf(frequency_array,tf_fbg2.transmission,'si',save_dispersion);
% Extract dispersion from transfer function.
% Calculate group delay and dispersion.
% We use the transfer function in transmission to avoid phase
% discontinuities. OK for symmetric apodisation profiles.


%--------------------------------------------------------------------------
% FBG3
%--------------------------------------------------------------------------
tf_fbg3 = opt_tf_fbg(freq,params_fbg3,numparams_fbg3);
% Calculate grating transfer function.

save_dispersion.status = 0;
save_dispersion.file_name = 'filter_tf.dat';
[phase_fbg3,freq_group_delay_fbg3,group_delay_fbg3,freq_dispersion_fbg3,dispersion_fbg3] = extract_dispersion_from_tf(frequency_array,tf_fbg3.transmission,'si',save_dispersion);
% Extract dispersion from transfer function.
% Calculate group delay and dispersion.
% We use the transfer function in transmission to avoid phase
% discontinuities. OK for symmetric apodisation profiles.


%--------------------------------------------------------------------------
% FBG4
%--------------------------------------------------------------------------
tf_fbg4 = opt_tf_fbg(freq,params_fbg4,numparams_fbg4);
% Calculate grating transfer function.

save_dispersion.status = 0;
save_dispersion.file_name = 'filter_tf.dat';
[phase_fbg4,freq_group_delay_fbg4,group_delay_fbg4,freq_dispersion_fbg4,dispersion_fbg4] = extract_dispersion_from_tf(frequency_array,tf_fbg4.transmission,'si',save_dispersion);
% Extract dispersion from transfer function.
% Calculate group delay and dispersion.
% We use the transfer function in transmission to avoid phase
% discontinuities. OK for symmetric apodisation profiles.



%--------------------------------------------------------------------------
% FBG5
%--------------------------------------------------------------------------
tf_fbg5 = opt_tf_fbg(freq,params_fbg5,numparams_fbg5);
% Calculate grating transfer function.

save_dispersion.status = 0;
save_dispersion.file_name = 'filter_tf.dat';
[phase_fbg5,freq_group_delay_fbg5,group_delay_fbg5,freq_dispersion_fbg5,dispersion_fbg5] = extract_dispersion_from_tf(frequency_array,tf_fbg5.reflection,'si',save_dispersion);
% Extract dispersion from transfer function.
% Calculate group delay and dispersion.


%--------------------------------------------------------------------------
% FBG6
%--------------------------------------------------------------------------
tf_fbg6 = opt_tf_fbg(freq,params_fbg6,numparams_fbg6);
% Calculate grating transfer function.

save_dispersion.status = 0;
save_dispersion.file_name = 'filter_tf.dat';
[phase_fbg6,freq_group_delay_fbg6,group_delay_fbg6,freq_dispersion_fbg6,dispersion_fbg6] = extract_dispersion_from_tf(frequency_array,tf_fbg6.reflection,'si',save_dispersion);
% Extract dispersion from transfer function.
% Calculate group delay and dispersion.





% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% PLOT FBG1
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
fig_name = 'tf_fbg1';
hfig = figure('Name',fig_name);

% -------------------------------------------------------------------------
% PLOT
% -------------------------------------------------------------------------
yyaxis left
h1 = plot(CONSTANT.c./freq/1.0e-9,abs(tf_fbg1.reflection).^2,'LineWidth',1.5,'Color','b','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on;
yyaxis right
h2 = plot(CONSTANT.c./(reference_frequency + freq_group_delay_fbg1)/1.0e-9,group_delay_fbg1/1.0e-12,'LineWidth',1.5,'Color','r','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');

% -------------------------------------------------------------------------
% LEGEND
% -------------------------------------------------------------------------
hleg = legend([h1,h2],{'',''});

hleg.Title.String = '';
hleg.Location = 'SouthEast';
hleg.Orientation = 'vertical';
hleg.NumColumns = 1;

hleg.Box = 'off';
hleg.EdgeColor = 'black';
hleg.TextColor = 'red';
hleg.Color = 'white';
hleg.LineWidth = 1;

hleg.Interpreter = fig.interpreter;
hleg.FontName = fig.font_name;
hleg.FontSize = fig.font_size;

hleg.Visible = 'on';


% -------------------------------------------------------------------------
% AXES
% -------------------------------------------------------------------------
ax = gca;

ax.LineWidth = 1.5;
% ax.XMinorTick = 'on';
% ax.YMinorTick = 'on';
ax.TickLength = [0.01 0.025]*1.5;
ax.TickLabelInterpreter = fig.interpreter;

ax.XAxis.FontName = fig.font_name;
ax.XAxis.FontSize = fig.font_size;
ax.XAxis.Color = 'k';
xlim([1549.6 1550.4]);
ax.XTick = [1549.6:0.2:1550.4];
% xtickformat(ax,'%3.3f')
% ax.XTickLabel = {'','',''}

% yyaxis left
ax.YAxis(1).FontName = fig.font_name;
ax.YAxis(1).FontSize = fig.font_size;
ax.YAxis(1).Color = 'k';
ylim([0 1.1]);
ax.YTick = [0:0.2:1];
% ytickformat(ax,'%3.3f')
% ax.YTickLabel = {'','',''}

yyaxis right
ax.YAxis(2).FontName = fig.font_name;
ax.YAxis(2).FontSize = fig.font_size;
ax.YAxis(2).Color = 'k';
ylim([0 100]);
ax.YTick = [0:20:100];
% ytickformat(ax,'%3.3f')
% ax.YTickLabel = {'','',''}







% -------------------------------------------------------------------------
% AXES LABELS
% -------------------------------------------------------------------------
x1 = xlabel('wavelength (nm)');
x1.Interpreter = fig.interpreter;
x1.FontName = fig.font_name;
x1.FontSize = fig.font_size;
x1.FontWeight = 'normal';
x1.Color = 'k';

yyaxis left
yl = ylabel('reflectivity');
yl.Interpreter = fig.interpreter;
yl.FontName = fig.font_name;
yl.FontSize = fig.font_size;
yl.FontWeight = 'normal';
yl.Color = 'k';

yyaxis right
y2 = ylabel('group delay (ps)');
y2.Interpreter = fig.interpreter;
y2.FontName = fig.font_name;
y2.FontSize = fig.font_size;
% y2.FontWeight = 'normal';
% y2.Color = 'k';



% -------------------------------------------------------------------------
% GRID LINES
% -------------------------------------------------------------------------
ax.XGrid = 0;
ax.XMinorGrid = 0;

ax.YGrid = 0;
ax.YMinorGrid = 0;

ax.GridColorMode = 'manual';
ax.GridColor = 'k';

ax.MinorGridColorMode = 'manual';
ax.MinorGridColor = 'k';

ax.GridAlphaMode = 'manual';
ax.GridAlpha = 0.25;

ax.MinorGridAlphaMode = 'manual';
ax.MinorGridAlpha = 0.25;

ax.GridLineStyle = '-';
ax.MinorGridLineStyle = '--';


% -------------------------------------------------------------------------
% ADD TEXT
% -------------------------------------------------------------------------
% tt = text(0.63,-0.5,'$x \simeq 0.6$');
% tt.Interpreter = fig.interpreter;
% tt.FontSize = fig.font_size;
% tt.FontName = fig.font_name;
% tt.FontWeight = 'normal';
% tt.Color = 'black';
% tt.HorizontalAlignment = 'left';


%--------------------------------------------------------------------------
% HORIZONTAL OR VERTICAL LINES
%--------------------------------------------------------------------------
% xline = 7.5;
% yy = get(gca,'ylim');
% plot([xline xline],yy,'k--','LineWidth',1);
% % Vertical line
% 
% yline = 0.4;
% xx = get(gca,'xlim');
% plot(xx,[yline yline],'k','Linewidth',1);
% % Horizontal line




% Not sure what to do with this so far...
% fig.PaperUnits = 'centimeters';
% fig.PaperType = 'a4';
% fig.PaperOrientation = 'landscape';
% fig.PaperPosition = [5 5 17 8.5]


% -------------------------------------------------------------------------
% PRINT TO FILE
% -------------------------------------------------------------------------
% print(fig_name,'-dmeta');

% print(fig_name,'-dpdf');
% Command =['pdfcrop ' fig_name '.pdf ' fig_name '.pdf'];
% system(Command);




% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% PLOT FBG2
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
fig_name = 'tf_fbg2';
hfig = figure('Name',fig_name);

% -------------------------------------------------------------------------
% PLOT
% -------------------------------------------------------------------------
yyaxis left
h1 = plot(CONSTANT.c./freq/1.0e-9,abs(tf_fbg2.reflection).^2,'LineWidth',1.5,'Color','b','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on;
yyaxis right
h2 = plot(CONSTANT.c./(reference_frequency + freq_group_delay_fbg2)/1.0e-9,group_delay_fbg2/1.0e-12,'LineWidth',1.5,'Color','r','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');

% -------------------------------------------------------------------------
% LEGEND
% -------------------------------------------------------------------------
hleg = legend([h1,h2],{'',''});

hleg.Title.String = '';
hleg.Location = 'SouthEast';
hleg.Orientation = 'vertical';
hleg.NumColumns = 1;

hleg.Box = 'off';
hleg.EdgeColor = 'black';
hleg.TextColor = 'red';
hleg.Color = 'white';
hleg.LineWidth = 1;

hleg.Interpreter = fig.interpreter;
hleg.FontName = fig.font_name;
hleg.FontSize = fig.font_size;

hleg.Visible = 'on';


% -------------------------------------------------------------------------
% AXES
% -------------------------------------------------------------------------
ax = gca;

ax.LineWidth = 1.5;
% ax.XMinorTick = 'on';
% ax.YMinorTick = 'on';
ax.TickLength = [0.01 0.025]*1.5;
ax.TickLabelInterpreter = fig.interpreter;

ax.XAxis.FontName = fig.font_name;
ax.XAxis.FontSize = fig.font_size;
ax.XAxis.Color = 'k';
xlim([1549.6 1550.4]);
ax.XTick = [1549.6:0.2:1550.4];
% xtickformat(ax,'%3.3f')
% ax.XTickLabel = {'','',''}

% yyaxis left
ax.YAxis(1).FontName = fig.font_name;
ax.YAxis(1).FontSize = fig.font_size;
ax.YAxis(1).Color = 'k';
ylim([0 1.1]);
ax.YTick = [0:0.2:1];
% ytickformat(ax,'%3.3f')
% ax.YTickLabel = {'','',''}

yyaxis right
ax.YAxis(2).FontName = fig.font_name;
ax.YAxis(2).FontSize = fig.font_size;
ax.YAxis(2).Color = 'k';
ylim([-100 500]);
ax.YTick = [-100:100:500];
% ytickformat(ax,'%3.3f')
% ax.YTickLabel = {'','',''}


% -------------------------------------------------------------------------
% AXES LABELS
% -------------------------------------------------------------------------
x1 = xlabel('wavelength (nm)');
x1.Interpreter = fig.interpreter;
x1.FontName = fig.font_name;
x1.FontSize = fig.font_size;
x1.FontWeight = 'normal';
x1.Color = 'k';

yyaxis left
yl = ylabel('reflectivity');
yl.Interpreter = fig.interpreter;
yl.FontName = fig.font_name;
yl.FontSize = fig.font_size;
yl.FontWeight = 'normal';
yl.Color = 'k';

yyaxis right
y2 = ylabel('group delay (ps)');
y2.Interpreter = fig.interpreter;
y2.FontName = fig.font_name;
y2.FontSize = fig.font_size;
% y2.FontWeight = 'normal';
% y2.Color = 'k';



% -------------------------------------------------------------------------
% GRID LINES
% -------------------------------------------------------------------------
ax.XGrid = 0;
ax.XMinorGrid = 0;

ax.YGrid = 0;
ax.YMinorGrid = 0;

ax.GridColorMode = 'manual';
ax.GridColor = 'k';

ax.MinorGridColorMode = 'manual';
ax.MinorGridColor = 'k';

ax.GridAlphaMode = 'manual';
ax.GridAlpha = 0.25;

ax.MinorGridAlphaMode = 'manual';
ax.MinorGridAlpha = 0.25;

ax.GridLineStyle = '-';
ax.MinorGridLineStyle = '--';


% -------------------------------------------------------------------------
% ADD TEXT
% -------------------------------------------------------------------------
% tt = text(0.63,-0.5,'$x \simeq 0.6$');
% tt.Interpreter = fig.interpreter;
% tt.FontSize = fig.font_size;
% tt.FontName = fig.font_name;
% tt.FontWeight = 'normal';
% tt.Color = 'black';
% tt.HorizontalAlignment = 'left';


%--------------------------------------------------------------------------
% HORIZONTAL OR VERTICAL LINES
%--------------------------------------------------------------------------
% xline = 7.5;
% yy = get(gca,'ylim');
% plot([xline xline],yy,'k--','LineWidth',1);
% % Vertical line
% 
% yline = 0.4;
% xx = get(gca,'xlim');
% plot(xx,[yline yline],'k','Linewidth',1);
% % Horizontal line




% Not sure what to do with this so far...
% fig.PaperUnits = 'centimeters';
% fig.PaperType = 'a4';
% fig.PaperOrientation = 'landscape';
% fig.PaperPosition = [5 5 17 8.5]


% -------------------------------------------------------------------------
% PRINT TO FILE
% -------------------------------------------------------------------------
% print(fig_name,'-dmeta');

% print(fig_name,'-dpdf');
% Command =['pdfcrop ' fig_name '.pdf ' fig_name '.pdf'];
% system(Command);









% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% PLOT FBG3
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
fig_name = 'tf_fbg3';
hfig = figure('Name',fig_name);

% -------------------------------------------------------------------------
% PLOT
% -------------------------------------------------------------------------
yyaxis left
h1 = plot(CONSTANT.c./freq/1.0e-9,abs(tf_fbg3.reflection).^2,'LineWidth',1.5,'Color','b','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on;
yyaxis right
h2 = plot(CONSTANT.c./(reference_frequency + freq_group_delay_fbg3)/1.0e-9,group_delay_fbg3/1.0e-12,'LineWidth',1.5,'Color','r','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');

% -------------------------------------------------------------------------
% LEGEND
% -------------------------------------------------------------------------
hleg = legend([h1,h2],{'',''});

hleg.Title.String = '';
hleg.Location = 'SouthEast';
hleg.Orientation = 'vertical';
hleg.NumColumns = 1;

hleg.Box = 'off';
hleg.EdgeColor = 'black';
hleg.TextColor = 'red';
hleg.Color = 'white';
hleg.LineWidth = 1;

hleg.Interpreter = fig.interpreter;
hleg.FontName = fig.font_name;
hleg.FontSize = fig.font_size;

hleg.Visible = 'on';


% -------------------------------------------------------------------------
% AXES
% -------------------------------------------------------------------------
ax = gca;

ax.LineWidth = 1.5;
% ax.XMinorTick = 'on';
% ax.YMinorTick = 'on';
ax.TickLength = [0.01 0.025]*1.5;
ax.TickLabelInterpreter = fig.interpreter;

ax.XAxis.FontName = fig.font_name;
ax.XAxis.FontSize = fig.font_size;
ax.XAxis.Color = 'k';
xlim([1549.6 1550.4]);
ax.XTick = [1549.6:0.2:1550.4];
% xtickformat(ax,'%3.3f')
% ax.XTickLabel = {'','',''}

% yyaxis left
ax.YAxis(1).FontName = fig.font_name;
ax.YAxis(1).FontSize = fig.font_size;
ax.YAxis(1).Color = 'k';
ylim([0 1.1]);
ax.YTick = [0:0.2:1];
% ytickformat(ax,'%3.3f')
% ax.YTickLabel = {'','',''}

yyaxis right
ax.YAxis(2).FontName = fig.font_name;
ax.YAxis(2).FontSize = fig.font_size;
ax.YAxis(2).Color = 'k';
ylim([-100 500]);
ax.YTick = [-100:100:500];
% ytickformat(ax,'%3.3f')
% ax.YTickLabel = {'','',''}

% -------------------------------------------------------------------------
% AXES LABELS
% -------------------------------------------------------------------------
x1 = xlabel('wavelength (nm)');
x1.Interpreter = fig.interpreter;
x1.FontName = fig.font_name;
x1.FontSize = fig.font_size;
x1.FontWeight = 'normal';
x1.Color = 'k';

yyaxis left
yl = ylabel('reflectivity');
yl.Interpreter = fig.interpreter;
yl.FontName = fig.font_name;
yl.FontSize = fig.font_size;
yl.FontWeight = 'normal';
yl.Color = 'k';

yyaxis right
y2 = ylabel('group delay (ps)');
y2.Interpreter = fig.interpreter;
y2.FontName = fig.font_name;
y2.FontSize = fig.font_size;
% y2.FontWeight = 'normal';
% y2.Color = 'k';



% -------------------------------------------------------------------------
% GRID LINES
% -------------------------------------------------------------------------
ax.XGrid = 0;
ax.XMinorGrid = 0;

ax.YGrid = 0;
ax.YMinorGrid = 0;

ax.GridColorMode = 'manual';
ax.GridColor = 'k';

ax.MinorGridColorMode = 'manual';
ax.MinorGridColor = 'k';

ax.GridAlphaMode = 'manual';
ax.GridAlpha = 0.25;

ax.MinorGridAlphaMode = 'manual';
ax.MinorGridAlpha = 0.25;

ax.GridLineStyle = '-';
ax.MinorGridLineStyle = '--';


% -------------------------------------------------------------------------
% ADD TEXT
% -------------------------------------------------------------------------
% tt = text(0.63,-0.5,'$x \simeq 0.6$');
% tt.Interpreter = fig.interpreter;
% tt.FontSize = fig.font_size;
% tt.FontName = fig.font_name;
% tt.FontWeight = 'normal';
% tt.Color = 'black';
% tt.HorizontalAlignment = 'left';


%--------------------------------------------------------------------------
% HORIZONTAL OR VERTICAL LINES
%--------------------------------------------------------------------------
% xline = 7.5;
% yy = get(gca,'ylim');
% plot([xline xline],yy,'k--','LineWidth',1);
% % Vertical line
% 
% yline = 0.4;
% xx = get(gca,'xlim');
% plot(xx,[yline yline],'k','Linewidth',1);
% % Horizontal line




% Not sure what to do with this so far...
% fig.PaperUnits = 'centimeters';
% fig.PaperType = 'a4';
% fig.PaperOrientation = 'landscape';
% fig.PaperPosition = [5 5 17 8.5]


% -------------------------------------------------------------------------
% PRINT TO FILE
% -------------------------------------------------------------------------
% print(fig_name,'-dmeta');

% print(fig_name,'-dpdf');
% Command =['pdfcrop ' fig_name '.pdf ' fig_name '.pdf'];
% system(Command);





% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% PLOT FBG4
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
fig_name = 'tf_fbg4';
hfig = figure('Name',fig_name);

% -------------------------------------------------------------------------
% PLOT
% -------------------------------------------------------------------------
yyaxis left
h1 = plot(CONSTANT.c./freq/1.0e-9,abs(tf_fbg4.reflection).^2,'LineWidth',1.5,'Color','b','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on;
yyaxis right
h2 = plot(CONSTANT.c./(reference_frequency + freq_group_delay_fbg4)/1.0e-9,group_delay_fbg4/1.0e-12,'LineWidth',1.5,'Color','r','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');

% -------------------------------------------------------------------------
% LEGEND
% -------------------------------------------------------------------------
hleg = legend([h1,h2],{'',''});

hleg.Title.String = '';
hleg.Location = 'SouthEast';
hleg.Orientation = 'vertical';
hleg.NumColumns = 1;

hleg.Box = 'off';
hleg.EdgeColor = 'black';
hleg.TextColor = 'red';
hleg.Color = 'white';
hleg.LineWidth = 1;

hleg.Interpreter = fig.interpreter;
hleg.FontName = fig.font_name;
hleg.FontSize = fig.font_size;

hleg.Visible = 'on';


% -------------------------------------------------------------------------
% AXES
% -------------------------------------------------------------------------
ax = gca;

ax.LineWidth = 1.5;
% ax.XMinorTick = 'on';
% ax.YMinorTick = 'on';
ax.TickLength = [0.01 0.025]*1.5;
ax.TickLabelInterpreter = fig.interpreter;

ax.XAxis.FontName = fig.font_name;
ax.XAxis.FontSize = fig.font_size;
ax.XAxis.Color = 'k';
xlim([1549.6 1550.4]);
ax.XTick = [1549.6:0.2:1550.4];
% xtickformat(ax,'%3.3f')
% ax.XTickLabel = {'','',''}

% yyaxis left
ax.YAxis(1).FontName = fig.font_name;
ax.YAxis(1).FontSize = fig.font_size;
ax.YAxis(1).Color = 'k';
ylim([0 1.1]);
ax.YTick = [0:0.2:1];
% ytickformat(ax,'%3.3f')
% ax.YTickLabel = {'','',''}

yyaxis right
ax.YAxis(2).FontName = fig.font_name;
ax.YAxis(2).FontSize = fig.font_size;
ax.YAxis(2).Color = 'k';
ylim([-100 500]);
ax.YTick = [-100:100:500];
% ytickformat(ax,'%3.3f')
% ax.YTickLabel = {'','',''}


% -------------------------------------------------------------------------
% AXES LABELS
% -------------------------------------------------------------------------
x1 = xlabel('wavelength (nm)');
x1.Interpreter = fig.interpreter;
x1.FontName = fig.font_name;
x1.FontSize = fig.font_size;
x1.FontWeight = 'normal';
x1.Color = 'k';

yyaxis left
yl = ylabel('reflectivity');
yl.Interpreter = fig.interpreter;
yl.FontName = fig.font_name;
yl.FontSize = fig.font_size;
yl.FontWeight = 'normal';
yl.Color = 'k';

yyaxis right
y2 = ylabel('group delay (ps)');
y2.Interpreter = fig.interpreter;
y2.FontName = fig.font_name;
y2.FontSize = fig.font_size;
% y2.FontWeight = 'normal';
% y2.Color = 'k';



% -------------------------------------------------------------------------
% GRID LINES
% -------------------------------------------------------------------------
ax.XGrid = 0;
ax.XMinorGrid = 0;

ax.YGrid = 0;
ax.YMinorGrid = 0;

ax.GridColorMode = 'manual';
ax.GridColor = 'k';

ax.MinorGridColorMode = 'manual';
ax.MinorGridColor = 'k';

ax.GridAlphaMode = 'manual';
ax.GridAlpha = 0.25;

ax.MinorGridAlphaMode = 'manual';
ax.MinorGridAlpha = 0.25;

ax.GridLineStyle = '-';
ax.MinorGridLineStyle = '--';


% -------------------------------------------------------------------------
% ADD TEXT
% -------------------------------------------------------------------------
% tt = text(0.63,-0.5,'$x \simeq 0.6$');
% tt.Interpreter = fig.interpreter;
% tt.FontSize = fig.font_size;
% tt.FontName = fig.font_name;
% tt.FontWeight = 'normal';
% tt.Color = 'black';
% tt.HorizontalAlignment = 'left';


%--------------------------------------------------------------------------
% HORIZONTAL OR VERTICAL LINES
%--------------------------------------------------------------------------
% xline = 7.5;
% yy = get(gca,'ylim');
% plot([xline xline],yy,'k--','LineWidth',1);
% % Vertical line
% 
% yline = 0.4;
% xx = get(gca,'xlim');
% plot(xx,[yline yline],'k','Linewidth',1);
% % Horizontal line




% Not sure what to do with this so far...
% fig.PaperUnits = 'centimeters';
% fig.PaperType = 'a4';
% fig.PaperOrientation = 'landscape';
% fig.PaperPosition = [5 5 17 8.5]


% -------------------------------------------------------------------------
% PRINT TO FILE
% -------------------------------------------------------------------------
% print(fig_name,'-dmeta');

% print(fig_name,'-dpdf');
% Command =['pdfcrop ' fig_name '.pdf ' fig_name '.pdf'];
% system(Command);






% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% PLOT FBG5
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
fig_name = 'tf_fbg5';
hfig = figure('Name',fig_name);

% -------------------------------------------------------------------------
% PLOT
% -------------------------------------------------------------------------
yyaxis left
h1 = plot(CONSTANT.c./freq/1.0e-9,abs(tf_fbg5.reflection).^2,'LineWidth',1.5,'Color','b','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on;
yyaxis right
h2 = plot(CONSTANT.c./(reference_frequency + freq_group_delay_fbg5)/1.0e-9,group_delay_fbg5/1.0e-12,'LineWidth',1.5,'Color','r','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');

% -------------------------------------------------------------------------
% LEGEND
% -------------------------------------------------------------------------
hleg = legend([h1,h2],{'',''});

hleg.Title.String = '';
hleg.Location = 'SouthEast';
hleg.Orientation = 'vertical';
hleg.NumColumns = 1;

hleg.Box = 'off';
hleg.EdgeColor = 'black';
hleg.TextColor = 'red';
hleg.Color = 'white';
hleg.LineWidth = 1;

hleg.Interpreter = fig.interpreter;
hleg.FontName = fig.font_name;
hleg.FontSize = fig.font_size;

hleg.Visible = 'on';


% -------------------------------------------------------------------------
% AXES
% -------------------------------------------------------------------------
ax = gca;

ax.LineWidth = 1.5;
% ax.XMinorTick = 'on';
% ax.YMinorTick = 'on';
ax.TickLength = [0.01 0.025]*1.5;
ax.TickLabelInterpreter = fig.interpreter;

ax.XAxis.FontName = fig.font_name;
ax.XAxis.FontSize = fig.font_size;
ax.XAxis.Color = 'k';
xlim([1548 1552]);
ax.XTick = [1548:1:1552];
% xtickformat(ax,'%3.3f')
% ax.XTickLabel = {'','',''}

% yyaxis left
ax.YAxis(1).FontName = fig.font_name;
ax.YAxis(1).FontSize = fig.font_size;
ax.YAxis(1).Color = 'k';
ylim([0 1.1]);
ax.YTick = [0:0.2:1];
% ytickformat(ax,'%3.3f')
% ax.YTickLabel = {'','',''}

yyaxis right
ax.YAxis(2).FontName = fig.font_name;
ax.YAxis(2).FontSize = fig.font_size;
ax.YAxis(2).Color = 'k';
ylim([-500 1500]);
ax.YTick = [-500:100:1500];
% ytickformat(ax,'%3.3f')
% ax.YTickLabel = {'','',''}


% -------------------------------------------------------------------------
% AXES LABELS
% -------------------------------------------------------------------------
x1 = xlabel('wavelength (nm)');
x1.Interpreter = fig.interpreter;
x1.FontName = fig.font_name;
x1.FontSize = fig.font_size;
x1.FontWeight = 'normal';
x1.Color = 'k';

yyaxis left
yl = ylabel('reflectivity');
yl.Interpreter = fig.interpreter;
yl.FontName = fig.font_name;
yl.FontSize = fig.font_size;
yl.FontWeight = 'normal';
yl.Color = 'k';

yyaxis right
y2 = ylabel('group delay (ps)');
y2.Interpreter = fig.interpreter;
y2.FontName = fig.font_name;
y2.FontSize = fig.font_size;
% y2.FontWeight = 'normal';
% y2.Color = 'k';



% -------------------------------------------------------------------------
% GRID LINES
% -------------------------------------------------------------------------
ax.XGrid = 0;
ax.XMinorGrid = 0;

ax.YGrid = 0;
ax.YMinorGrid = 0;

ax.GridColorMode = 'manual';
ax.GridColor = 'k';

ax.MinorGridColorMode = 'manual';
ax.MinorGridColor = 'k';

ax.GridAlphaMode = 'manual';
ax.GridAlpha = 0.25;

ax.MinorGridAlphaMode = 'manual';
ax.MinorGridAlpha = 0.25;

ax.GridLineStyle = '-';
ax.MinorGridLineStyle = '--';


% -------------------------------------------------------------------------
% ADD TEXT
% -------------------------------------------------------------------------
% tt = text(0.63,-0.5,'$x \simeq 0.6$');
% tt.Interpreter = fig.interpreter;
% tt.FontSize = fig.font_size;
% tt.FontName = fig.font_name;
% tt.FontWeight = 'normal';
% tt.Color = 'black';
% tt.HorizontalAlignment = 'left';


%--------------------------------------------------------------------------
% HORIZONTAL OR VERTICAL LINES
%--------------------------------------------------------------------------
% xline = 7.5;
% yy = get(gca,'ylim');
% plot([xline xline],yy,'k--','LineWidth',1);
% % Vertical line
% 
% yline = 0.4;
% xx = get(gca,'xlim');
% plot(xx,[yline yline],'k','Linewidth',1);
% % Horizontal line




% Not sure what to do with this so far...
% fig.PaperUnits = 'centimeters';
% fig.PaperType = 'a4';
% fig.PaperOrientation = 'landscape';
% fig.PaperPosition = [5 5 17 8.5]


% -------------------------------------------------------------------------
% PRINT TO FILE
% -------------------------------------------------------------------------
% print(fig_name,'-dmeta');

% print(fig_name,'-dpdf');
% Command =['pdfcrop ' fig_name '.pdf ' fig_name '.pdf'];
% system(Command);




% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% PLOT FBG6
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
fig_name = 'tf_fbg6';
hfig = figure('Name',fig_name);

% -------------------------------------------------------------------------
% PLOT
% -------------------------------------------------------------------------
yyaxis left
h1 = plot(CONSTANT.c./freq/1.0e-9,abs(tf_fbg6.reflection).^2,'LineWidth',1.5,'Color','b','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on;
yyaxis right
h2 = plot(CONSTANT.c./(reference_frequency + freq_group_delay_fbg6)/1.0e-9,group_delay_fbg6/1.0e-12,'LineWidth',1.5,'Color','r','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');

% -------------------------------------------------------------------------
% LEGEND
% -------------------------------------------------------------------------
hleg = legend([h1,h2],{'',''});

hleg.Title.String = '';
hleg.Location = 'SouthEast';
hleg.Orientation = 'vertical';
hleg.NumColumns = 1;

hleg.Box = 'off';
hleg.EdgeColor = 'black';
hleg.TextColor = 'red';
hleg.Color = 'white';
hleg.LineWidth = 1;

hleg.Interpreter = fig.interpreter;
hleg.FontName = fig.font_name;
hleg.FontSize = fig.font_size;

hleg.Visible = 0;


% -------------------------------------------------------------------------
% AXES
% -------------------------------------------------------------------------
ax = gca;

ax.LineWidth = 1.5;
% ax.XMinorTick = 'on';
% ax.YMinorTick = 'on';
ax.TickLength = [0.01 0.025]*1.5;
ax.TickLabelInterpreter = fig.interpreter;

ax.XAxis.FontName = fig.font_name;
ax.XAxis.FontSize = fig.font_size;
ax.XAxis.Color = 'k';
xlim([1548 1552]);
ax.XTick = [1548:1:1552];
% xtickformat(ax,'%3.3f')
% ax.XTickLabel = {'','',''}

% yyaxis left
ax.YAxis(1).FontName = fig.font_name;
ax.YAxis(1).FontSize = fig.font_size;
ax.YAxis(1).Color = 'k';
ylim([0 1.1]);
ax.YTick = [0:0.2:1];
% ytickformat(ax,'%3.3f')
% ax.YTickLabel = {'','',''}

yyaxis right
ax.YAxis(2).FontName = fig.font_name;
ax.YAxis(2).FontSize = fig.font_size;
ax.YAxis(2).Color = 'k';
ylim([-500 1500]);
ax.YTick = [-500:500:1500];
% ytickformat(ax,'%3.3f')
% ax.YTickLabel = {'','',''}


% -------------------------------------------------------------------------
% AXES LABELS
% -------------------------------------------------------------------------
x1 = xlabel('wavelength (nm)');
x1.Interpreter = fig.interpreter;
x1.FontName = fig.font_name;
x1.FontSize = fig.font_size;
x1.FontWeight = 'normal';
x1.Color = 'k';

yyaxis left
yl = ylabel('reflectivity');
yl.Interpreter = fig.interpreter;
yl.FontName = fig.font_name;
yl.FontSize = fig.font_size;
yl.FontWeight = 'normal';
yl.Color = 'k';

yyaxis right
y2 = ylabel('group delay (ps)');
y2.Interpreter = fig.interpreter;
y2.FontName = fig.font_name;
y2.FontSize = fig.font_size;
% y2.FontWeight = 'normal';
% y2.Color = 'k';



% -------------------------------------------------------------------------
% GRID LINES
% -------------------------------------------------------------------------
ax.XGrid = 0;
ax.XMinorGrid = 0;

ax.YGrid = 0;
ax.YMinorGrid = 0;

ax.GridColorMode = 'manual';
ax.GridColor = 'k';

ax.MinorGridColorMode = 'manual';
ax.MinorGridColor = 'k';

ax.GridAlphaMode = 'manual';
ax.GridAlpha = 0.25;

ax.MinorGridAlphaMode = 'manual';
ax.MinorGridAlpha = 0.25;

ax.GridLineStyle = '-';
ax.MinorGridLineStyle = '--';


% -------------------------------------------------------------------------
% ADD TEXT
% -------------------------------------------------------------------------
% tt = text(0.63,-0.5,'$x \simeq 0.6$');
% tt.Interpreter = fig.interpreter;
% tt.FontSize = fig.font_size;
% tt.FontName = fig.font_name;
% tt.FontWeight = 'normal';
% tt.Color = 'black';
% tt.HorizontalAlignment = 'left';


%--------------------------------------------------------------------------
% HORIZONTAL OR VERTICAL LINES
%--------------------------------------------------------------------------
% xline = 7.5;
% yy = get(gca,'ylim');
% plot([xline xline],yy,'k--','LineWidth',1);
% % Vertical line
% 
% yline = 0.4;
% xx = get(gca,'xlim');
% plot(xx,[yline yline],'k','Linewidth',1);
% % Horizontal line




% Not sure what to do with this so far...
% fig.PaperUnits = 'centimeters';
% fig.PaperType = 'a4';
% fig.PaperOrientation = 'landscape';
% fig.PaperPosition = [5 5 17 8.5]


% -------------------------------------------------------------------------
% PRINT TO FILE
% -------------------------------------------------------------------------
% print(fig_name,'-dmeta');

% print(fig_name,'-dpdf');
% Command =['pdfcrop ' fig_name '.pdf ' fig_name '.pdf'];
% system(Command);


% -------------------------------------------------------------------------
% Figure alignment
% -------------------------------------------------------------------------
% Requires alignfigs function
% Available at https://github.com/nickhale/alignfigs
% alignfigs(2)

%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Save data
% -------------------------------------------------------------------------
% ------------------------------------------------------------------------- 
% file_name_data = [file_name_core_data '_' time_stamp]
% save(file_name_data,'');
        
        


% ------------------------------------------------------------------------- 
% ------------------------------------------------------------------------- 
% Display duration
% ------------------------------------------------------------------------- 
% ------------------------------------------------------------------------- 
core_display_duration(start_time,datetime("now"));


