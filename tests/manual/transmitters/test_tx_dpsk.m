% -------------------------------------------------------------------------
% Test of tx_dpsk.m
%
% Generation of OOK signals with 33% RZ, 50% RZ, 67% RZ and NRZ pulse
% shapes.
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
% Sequences generation
% -------------------------------------------------------------------------
params_prbs.type = 'shift_register';%'de_bruijn';
% params_prbs.order = 7;
% params_prbs.poly = [7 6 0];%[7 3 0];[7 1 0];
% params_prbs.seed = [1 1 1 0 1 1 0];

params_prbs.order = 11;
params_prbs.poly = [11 9 0];%[11 8 5 2 0];
params_prbs.seed = [0 1 1 0 1 1 1 0 0 1 0];
bit_pattern = logical_prbs(params_prbs);
bit_pattern = logical_adapt_binary_sequence(bit_pattern,nsymbols);
% Generate random sequence and adapt length to simulation time window

bit_pattern_diff_enc = logical_differential_encoder_binary(bit_pattern);
% Differentially-encoded sequence


% -------------------------------------------------------------------------
% Electrical signal 
% -------------------------------------------------------------------------
rise_time = 1/symbol_rate/4;  
% NRZ signal rise time, in s

nrz_data_sig = elec_pulse_sequence_nrz(bit_pattern,rise_time);
% Encode the binary bit sequence onto an NRZ signal

nrz_data_sig_diff_enc = elec_pulse_sequence_nrz(bit_pattern_diff_enc,rise_time);
% Encode the differentially binary bit pattern onto an NRZ signal.


params_scope.visualisers = {'amplitude'};
params_scope.display_interval = [0 32/symbol_rate];
params_scope.save.emf = 0;
params_scope.name = 'Original NRZ data signal';
meas_scope(nrz_data_sig,params_scope);
% Oscilloscope

params_scope.name = 'Differentially-encoded signal';
meas_scope(nrz_data_sig_diff_enc,params_scope);
% Oscilloscope



% -------------------------------------------------------------------------
% PM DPSK
% -------------------------------------------------------------------------
params_tx.emission_frequency = reference_frequency;
params_tx.power = 1.0e-3;
params_tx.linewidth = 0;
params_tx.bit_pattern = bit_pattern_diff_enc;
params_tx.rise_time = 1/symbol_rate/4;
sig_pm = tx_dpsk_nrz_pm(params_tx);
% Generate DPSK signal using PM

params_scope.visualisers = {'power'};
params_scope.display_interval = [0 32/symbol_rate];
params_scope.save.emf = 0;
params_scope.name = 'Power - PM-DPSK';
meas_scope(sig_pm.x,params_scope);
% Oscilloscope

params_osa.pol = 'x';%'y','both';
params_osa.display_interval = [-6*symbol_rate 6*symbol_rate];
params_osa.resolution_bandwidth = 1e9;
params_osa.sensitivity = -60;
params_osa.display = 1;
params_osa.save.txt = 0;
params_osa.save.emf = 0;
params_osa.save.jpg = 0;
params_osa.name = 'Optical spectrum - PM-DPSK';
meas_osa(sig_pm,params_osa);
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
params_eye.name = 'Eye diagram - PM-DPSK';
meas_eye(sig_pm,params_eye);
% Eye diagram


params_constellation.type = 'scatter_transitions';%'transitions';%'scatter';
params_constellation.pol = 'x';% 'y','both';
params_constellation.nsamples_per_symbol = nsamples_per_symbol;
params_constellation.sampling_time = nsamples_per_symbol/2;
params_constellation.normalisation = 'mean';%'max';
params_constellation.save.emf = 0;
% params_constellation.line_color = 'b';
% params_constellation.line_width = 3;
% params_constellation.marker_type = 'o';
% params_constellation.marker_color = 'r';
% params_constellation.marker_size = 36;
% params_constellation.plot_circular_grid = 1;
% params_constellation.plot_axes = 0;
params_constellation.name = 'Constellation - PM-DPSK';
char_opt_constellation(sig_pm,params_constellation); 
% Plot optical constellation


%%
% -------------------------------------------------------------------------
% MZM DPSK
% -------------------------------------------------------------------------
sig_mzm = tx_dpsk_nrz_mzm(params_tx);
% Generate DPSK signal using MZM

params_scope.name = 'Power - MZM-DPSK';
meas_scope(sig_mzm.x,params_scope);
% Oscilloscope

params_osa.name = 'Optical spectrum - MZM-DPSK';
meas_osa(sig_mzm,params_osa);
% Optical spectrum analyser

params_eye.name = 'Eye diagram - MZM-DPSK';
meas_eye(sig_mzm,params_eye);
% Eye diagram

params_constellation.name = 'Constellation - MZM-DPSK';
char_opt_constellation(sig_mzm,params_constellation);
% Plot optical constellation

%%
% -------------------------------------------------------------------------
% Interferometric detection
% -------------------------------------------------------------------------
sig = sig_mzm;
sig = sig_pm;
% Decide which signal to detect

params_rx.type = 'id';
params_rx.obpf.type = 'gaussian';
params_rx.obpf.order = 1;
params_rx.obpf.bandwidth = 4*symbol_rate;
params_rx.obpf.centre_frequency = 0;
params_rx.mzdi.input_port = 'upper';
params_rx.mzdi.mode = 'tuned';
params_rx.mzdi.delay_target = 1/symbol_rate;
params_rx.mzdi.coupling_ratio_in = 0.5;
params_rx.mzdi.coupling_ratio_out = 0.5;
params_rx.mzdi.phase_shift = 0;
params_rx.upper.pd.responsivity = 1;
params_rx.upper.pd.shot_noise = 0;
params_rx.upper.pd.thermal_noise_density = 0;
params_rx.upper.pd.dark_current = 0;
params_rx.lower.pd.responsivity = 1;
params_rx.lower.pd.shot_noise = 0;
params_rx.lower.pd.thermal_noise_density = 0;
params_rx.lower.pd.dark_current = 0;
params_rx.upper.elpf.type = 'bessel';
params_rx.upper.elpf.order = 4;
params_rx.upper.elpf.f3dB = 0.75*symbol_rate;
params_rx.lower.elpf.type = 'bessel';
params_rx.lower.elpf.order = 4;
params_rx.lower.elpf.f3dB = 0.75*symbol_rate;
params_rx.bd.detection_mode = 'balanced';%'single_ended_upper';%'single_ended_lower';
params_rx.bd.electrical_delay = 0;
params_rx.bd.polarity = 1;%-1.

sig = rx_dd(sig,params_rx);

% For MZDI input_port = 'upper' and mode = 'tuned'
% 'single_ended_upper' -> compare to original data sequence and ignore 
% 1st bit. The optical signal is AMI.
% 'single_ended_lower' -> compare to inverted original data sequence and 
% ignore 1st bit. The optical signal is DB.


params_scope.visualisers = {'amplitude'};
params_scope.name = 'Detected signal';
meas_scope(sig,params_scope);
% Oscilloscope

params_eye.name = 'Detected eye diagram';
meas_eye(sig,params_eye);
% Eye diagram.

% We also check the tuning of the MZDI

sig_test = opt_broadband_time(0); 
% Generate white spectrum

params_osa.name = 'MZDI test - input spectrum';
meas_osa(sig_test,params_osa);

params_rx.mzdi.delay_target = 1/symbol_rate;
params_rx.mzdi.coupling_ratio_in = 0.5;
params_rx.mzdi.coupling_ratio_out = 0.5;
params_rx.mzdi.mode = 'tuned';
params_mzdi.phase_shift = 0;

[sig21,sig22] = opt_mzdi(sig_test,opt_nosig,params_rx.mzdi); 

params_osa.name = 'MZDI test - upper output port';
meas_osa(sig21,params_osa);

params_osa.name = 'MZDI test - lower output port';
meas_osa(sig22,params_osa);



% -------------------------------------------------------------------------
% Figure alignment
% -------------------------------------------------------------------------
% Requires alignfigs function
% Available at https://github.com/nickhale/alignfigs
alignfigs(3)

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