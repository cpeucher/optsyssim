% -------------------------------------------------------------------------
% Test of PAM4 signal generation
% 
% No DAC, analog summation of 2 electrical binary signals, which are then
% applied to a Mach-Zehnder modulator.
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
do_debug = 0;
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
nsamples_per_symbol = 64;
nsymbols = 2048;
symbol_rate = 10e9;

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
% Electrical signal 
% -------------------------------------------------------------------------
rise_time = 1/symbol_rate/4;  
% NRZ signal rise time, in s

params_prbs_1.type = 'shift_register';
params_prbs_1.order = 11;
params_prbs_1.poly = [11 9 0];%[11 8 5 2 0];
params_prbs_1.seed = [0 1 1 0 1 1 1 0 0 1 0];

params_prbs_2.type = 'shift_register';
params_prbs_2.order = 11;
params_prbs_2.poly = [11 8 5 2 0];%[11 9 0];
params_prbs_2.seed = [0 0 1 0 1 1 1 1 1 0 0];

bit_pattern_1 = logical_prbs(params_prbs_1);
bit_pattern_1 = logical_adapt_binary_sequence(bit_pattern_1,nsymbols);

bit_pattern_2 = logical_prbs(params_prbs_2);
bit_pattern_2 = logical_adapt_binary_sequence(bit_pattern_2,nsymbols);
% Generate random sequence and adapt length to simulation time window.

nrz_data_sig_1 = elec_pulse_sequence_nrz(bit_pattern_1,rise_time);
nrz_data_sig_2 = elec_pulse_sequence_nrz(bit_pattern_2,rise_time);
% Encode the binary bit patterns onto NRZ signals

pam4_data_sig = nrz_data_sig_1 + 0.5*nrz_data_sig_2;

params_scope.visualisers = {'amplitude'};
params_scope.display_interval = [0 32/symbol_rate];
params_scope.save.emf = 0;
params_scope.name = 'Waveform - Electrical PAM4 signal';
meas_scope(pam4_data_sig,params_scope); 
% Scope



%%
% -------------------------------------------------------------------------
% Generation of optical PAM4 signal
% -------------------------------------------------------------------------
params_tx.emission_frequency = reference_frequency;
params_tx.linewidth = 0;
params_tx.power = 1.0e-3;
params_tx.bit_pattern_1 = bit_pattern_1;
params_tx.bit_pattern_2 = bit_pattern_2;
params_tx.rise_time = 1/symbol_rate/4;
sig = tx_pam4_mzm(params_tx); 
% OOK RZ 33% transmitter

params_scope.visualisers = {'power'};
params_scope.display_interval = [0 32/symbol_rate];
params_scope.save.emf = 0;
params_scope.name = 'Waveform - PAM4';
meas_scope(sig.x,params_scope);
% Oscilloscope

params_osa.pol = 'x';%'y','both';
params_osa.display_interval = [-6*symbol_rate 6*symbol_rate];
params_osa.resolution_bandwidth = 1e9;
params_osa.sensitivity = -60;
params_osa.display = 1;
params_osa.save.txt = 0;
params_osa.save.emf = 0;
params_osa.save.jpg = 0;
params_osa.name = 'Optical spectrum - PAM4';
meas_osa(sig,params_osa);
% Optical spectrum analyser

params_eye.pol = 'x';%'y','both';
params_eye.neyes = 2;
params_eye.nsamples_per_symbol = nsamples_per_symbol;
params_eye.save.txt = 0;
params_eye.save.emf = 0;
params_eye.save.jpg = 0;
params_eye.save.vertical_scale_type = 'auto';%'fixed';
params_eye.save.vertical_scale_max = 1.0e-3;
params_eye.save.display_baseline = 1;
params_eye.colour_grade = 0;
params_eye.name = 'Eye diagram - PAM4';
meas_eye(sig,params_eye);
% Eye diagram




% -------------------------------------------------------------------------
% Figure alignment
% -------------------------------------------------------------------------
% Requires alignfigs function
% Available at https://github.com/nickhale/alignfigs
% alignfigs(3)

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