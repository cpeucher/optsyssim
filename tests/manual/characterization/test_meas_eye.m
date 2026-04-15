% -------------------------------------------------------------------------
% Test of meas_eye function 
% 
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
set(groot, "defaultFigurePosition", [680 458 560 420])
% Temporary fix so that the figures generated in R2025a have the same
% appearance as those generated in R2024b and earlier.
% See: https://se.mathworks.com/matlabcentral/answers/2175629-how-to-revert-the-figure-behavior-in-matlab-r2025a-and-newer-to-the-r2024b-style

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
nsymbols = 1024;
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

% -------------------------------------------------------------------------
% PRBS
% -------------------------------------------------------------------------
params_prbs.type = 'shift_register';%'de_bruijn';

% params_prbs.order = 7;
% params_prbs.poly = [7 6 0];%[7 3 0];[7 1 0];
% params_prbs.seed = [1 1 1 0 1 1 0];

params_prbs.order = 10;
params_prbs.poly = [10 3 0];%[10 7 0];
params_prbs.seed = [0 1 1 0 1 1 1 0 0 1];

% -------------------------------------------------------------------------
% Electrical signal 
% -------------------------------------------------------------------------
rise_time = 1/symbol_rate/4;  % NRZ signal rise time, in s.


        
        
% -------------------------------------------------------------------------        
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Here we go
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------    
% -------------------------------------------------------------------------

bit_pattern = logical_prbs(params_prbs);
bit_pattern = logical_adapt_binary_sequence(bit_pattern,nsymbols);
% Generate random sequence and adapt length to simulation time window


nrz_data_sig = elec_pulse_sequence_nrz(bit_pattern,rise_time);
% Encode the binary bit pattern onto an NRZ signal     

nbits_display = 20;
% Number of bits to display in waveform plots
      
params_scope.visualisers = {'amplitude'};
params_scope.display_interval = [0 nbits_display/symbol_rate];
params_scope.save.emf = 0;
params_scope.name = 'Electrical input signal';
meas_scope(nrz_data_sig,params_scope);
% Oscilloscope

params_tx.type = 'ook_nrz';
params_tx.emission_frequency = reference_frequency;
params_tx.power = 1.0e-3;
params_tx.linewidth = 0;
params_tx.bit_pattern = bit_pattern;
% params_tx.symbol_rate = symbol_rate;
params_tx.rise_time = 1/symbol_rate/4;
sig_tx = tx(params_tx);
% Optical transmitter

params_eye.pol = 'x';%'y','both';
params_eye.neyes = 2;
params_eye.nsamples_per_symbol = nsamples_per_symbol;
params_eye.display = 1;
params_eye.save.txt = 0;
params_eye.save.emf = 0;
params_eye.save.jpg = 0;
params_eye.save.vertical_scale_type = 'auto';%'fixed';
params_eye.save.vertical_scale_max = 1.0e-3;
params_eye.save.display_baseline = 0;
params_eye.colour_grade = 0;
params_eye.name = 'Optical eye diagram at tx output';
meas_eye(sig_tx,params_eye);
% Eye diagram



params_scope.visualisers = {'power'};
params_scope.display_interval = [0 nbits_display/symbol_rate];
params_scope.save.emf = 0;
params_scope.name = 'Optical signal at tx output';
meas_scope(sig_tx.x,params_scope);
% Oscilloscope

params_pin.pd.responsivity = 1;
params_pin.pd.thermal_noise_density = 10e-12;
params_pin.pd.shot_noise = 0;
params_pin.pd.dark_current = 0;
params_pin.elpf.type = 'bessel';%'none';'gaussian';'rc';'rectangular';
params_pin.elpf.order = 4;
params_pin.elpf.f3dB = 0.7*symbol_rate;
sigel_b2b = rx_pin(sig_tx,params_pin);
% PIN receiver 

params_eye.colour_grade = 1;
params_eye.name = 'Electrical eye diagram at rx output';
meas_eye(sigel_b2b,params_eye);
% Eye diagram



%%
% -------------------------------------------------------------------------
% We also compare the "legacy" plot of the eye diagram using the
% 'no_overlap' option with the 'overlap' option.
% For this purpose we introduce some inter-symbol interference by
% propagation over a dispersive medium.
% -------------------------------------------------------------------------

params_dispersion.dispersion = 17;           % Dispersion in ps/nm/km
params_dispersion.dispersion_slope = 0.058;  % Dispersion slope in ps/nm^2/km
params_dispersion.dispersion_curvature = 0;  % Dispersion curvature, in ps/nm^3/km
params_dispersion.dispersion_spec_frequency = reference_frequency;
z = 2e3; 
sig_disp = opt_dispersion(sig_tx,params_dispersion,z); 
% Linear and lossless dispersive element

sigel = rx_pin(sig_disp,params_pin);
% PIN receiver

sig_type =  'elec';
[sigel,imax] = rx_resynchronise(sigel,bit_pattern,sig_type);
% Resync signal, in case it is needed

% Plot eye diagram with the 'no_overlap' and 'overlap' options:

params_eye.colour_grade = 0;
params_eye.name = 'Electical eye diagram after transmission - no overlap';

tic
meas_eye(sigel,params_eye,'no_overlap');
toc 

params_eye.name = 'Electical eye diagram after transmission - overlap';

tic
meas_eye(sigel,params_eye,'overlap');
toc

% For such a high number of bits, we do not see any difference, but the 
% 'overlap' option is significantly slower. 
% Change the number of symbols to a much lower value to appreciate the
% difference between the 'no_overlap' and 'overlap' options.


%%
% -------------------------------------------------------------------------
% We also check the params_eye.display = 0 option to return the eye
% traces without plotting the eye diagram in the meas_eye function. 
% We then plot the eye diagram with improved control outside of the 
% function.
% -------------------------------------------------------------------------

params_eye.display = 0;
[timebase,traces] = meas_eye(sigel,params_eye,'no_overlap');

fig_name = [file_name_core_figure '_plot_outside'];
hfig = figure('Name',fig_name);
plot(timebase/1.0e-12,traces/1.0e-3,'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
xlabel('time (ps)','Interpreter',fig.interpreter)
ylabel('current / voltage (a.u.) ','Interpreter',fig.interpreter)
% legend('','Location','NorthEast','Box','off','Interpreter',fig.interpreter,'FontSize',fig.font_size);
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = fig.font_size;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
% xlim([0 18])
ylim([-0.1 1.2])

print_figure(hfig,do_print,do_add_figsize_to_filename,margin_figure)





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
% Display simulation duration
% ------------------------------------------------------------------------- 
% ------------------------------------------------------------------------- 
core_display_duration(start_time,datetime("now"));
