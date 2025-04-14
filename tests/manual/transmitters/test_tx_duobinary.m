% -------------------------------------------------------------------------
% Test of optical duobinary transmitters
%
% Generation of delay-and-add as well as low-pass filtered optical
% duobinary signals.
%
% We visualise the waveforms (optical power), eye diagrams, and optical
% power spectra.
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
nsamples_per_symbol = 32;
nsymbols = 128;
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
% Check waveforms
% -------------------------------------------------------------------------
params_prbs.type = 'shift_register';%'de_bruijn';
params_prbs.order = 7;
params_prbs.poly = [7 6 0];%[7 3 0];[7 1 0];
params_prbs.seed = [1 1 1 0 1 1 0];
bit_pattern = logical_prbs(params_prbs);
bit_pattern = logical_adapt_binary_sequence(bit_pattern,nsymbols);
% Generate random sequence and adapt length to simulation time window


bit_pattern_precoded = logical_differential_encoder_binary(not(bit_pattern));
% Pre-coded sequence

rise_time = 1/symbol_rate/4;
% NRZ signal rise time, in s

nrz_data_sig = elec_pulse_sequence_nrz(bit_pattern,rise_time);
% Encode the binary bit pattern onto an NRZ signal
nrz_data_sig_precoded = elec_pulse_sequence_nrz(bit_pattern_precoded,rise_time);
% Encode the precoded binary bit pattern onto an NRZ signal

duobinary_data_sig_da = 0.5*(nrz_data_sig_precoded + circshift(nrz_data_sig_precoded,[0 nsamples_per_symbol]));
% Delay-and-add duobinary encoder

elpf.type = 'bessel';
elpf.order = 4;
elpf.f3dB = 0.25*symbol_rate;
duobinary_data_sig_lp = elec_elpf(nrz_data_sig_precoded,elpf);
% Low pass filtering for duobinary encoding.

      
params_scope.visualisers = {'amplitude'};
params_scope.display_interval = [0 32/symbol_rate];
params_scope.save.txt = 0;
params_scope.save.emf = 0;
params_scope.name = 'Original NRZ data signal';
meas_scope(nrz_data_sig,params_scope);
% Oscilloscope

params_scope.name = 'Precoded NRZ data signal';
meas_scope(nrz_data_sig_precoded,params_scope);
% Oscilloscope

params_scope.name = 'Duobinary-encoded data signal (delay-and-add)';
meas_scope(duobinary_data_sig_da,params_scope);
% Oscilloscope

params_scope.name = 'Duobinary-encoded data signal (low-pass filter)';
meas_scope(duobinary_data_sig_lp,params_scope);
% Oscilloscope



%%
% -------------------------------------------------------------------------
% Delay-and-add optical duobinary transmitter
% -------------------------------------------------------------------------
params_tx.emission_frequency = reference_frequency;
params_tx.power = 1.0e-3;
params_tx.linewidth = 0;
params_tx.bit_pattern = bit_pattern;
params_tx.rise_time = 1/symbol_rate/4;
sig = tx_duobinary_nrz_delay_add(params_tx);
% Delay-and-add duobinary transmitter


params_scope.visualisers = {'power'};
params_scope.name = 'Waveform - Duobinary (delay-and-add)';
meas_scope(sig.x,params_scope);
% Oscilloscope

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
params_eye.name = 'Eye diagram - Duobinary (delay-and-add)';
meas_eye(sig,params_eye);
% Eye diagram

params_osa.pol = 'x';%'y','both';
params_osa.display_interval = [-6*symbol_rate 6*symbol_rate];
params_osa.resolution_bandwidth = 12.5e9/4;
params_osa.sensitivity = -60;
params_osa.display = 1;
params_osa.save.txt = 0;
params_osa.save.emf = 0;
params_osa.save.jpg = 0;
params_osa.name = 'Optical spectrum - Duobinary (delay-and-add)';
meas_osa(sig,params_osa);
% Optical spectrum analyser



%%
% -------------------------------------------------------------------------
% Low-pass filter optical duobinary transmitter
% -------------------------------------------------------------------------
sig = tx_duobinary_nrz_low_pass(params_tx);
% Low-pass filter duobinary transmitter


params_scope.name = 'Waveform - Duobinary (low-pass filter)';
meas_scope(sig.x,params_scope);
% Oscilloscope

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
params_eye.name = 'Eye diagram - Duobinary (low-pass filter)';
meas_eye(sig,params_eye);
% Eye diagram

params_osa.name = 'Optical spectrum - Duobinary (low-pass filter)';
meas_osa(sig,params_osa);
% Optical spectrum analyser



% -------------------------------------------------------------------------
% Figure alignment
% -------------------------------------------------------------------------
% Requires alignfigs function
% Available at https://github.com/nickhale/alignfigs
alignfigs(3)


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