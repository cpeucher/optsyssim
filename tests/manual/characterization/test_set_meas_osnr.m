% -------------------------------------------------------------------------
% Test of set_osnr and meas_osnr functions
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
nsamples_per_symbol = 64;
nsymbols = 1024;
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

params_prbs.type = 'shift_register';
params_prbs.order = 10;
params_prbs.poly = [10 3 0];%[10 7 0];
params_prbs.seed = [0 1 1 0 1 1 1 0 0 1];
bit_pattern = logical_prbs(params_prbs);
bit_pattern = logical_adapt_binary_sequence(bit_pattern,nsymbols);
% Generate random sequence and adapt length to simulation time window 

params_tx.type = 'ook_nrz';
params_tx.emission_frequency = reference_frequency - 50e9;
params_tx.power = 1.0e-3;
params_tx.linewidth = 0;
params_tx.bit_pattern = bit_pattern;
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
params_eye.save.display_baseline = 1;
params_eye.colour_grade = 0;

params_eye.name = ['Noise free optical eye diagram'];
meas_eye(sig_tx,params_eye);
% Check eye diagram of noise-free signal


params_osa.pol = 'both';%'x';%'y';
params_osa.display_interval = [frequency_array(1) frequency_array(end)];
params_osa.resolution_bandwidth = 1e9;
params_osa.sensitivity = -100;
params_osa.display = 1;
params_osa.save.txt = 0;
params_osa.save.emf = 0;
params_osa.save.jpg = 0;
params_osa.name = 'Check spectrum of noise-free signal';
meas_osa(sig_tx,params_osa); 
text(-100,params_osa.sensitivity + 5,['RB = ' num2str(params_osa.resolution_bandwidth/1.0e9) ' GHz'])
% Check optical spectrum of noise-free signal



%%
% Noise addition

ps = char_opt_average_power(sig_tx);
ps_dbm = 10*log10(ps/1.0e-3);
% Calculate the average power of the noise free signal

fprintf('\n\n%s%3.2f%s','Signal average power: ',ps_dbm,' dBm')



osnr_target_db = 20;
noise_bw = 12.5e9;
sig = set_osnr(sig_tx,osnr_target_db,noise_bw,2);
% Add noise

params_obpf.type = 'gaussian';%'none','rectangular_ideal','rectangular';
params_obpf.centre_frequency = 0;
params_obpf.bandwidth = 4*symbol_rate;
params_obpf.order = 4;
% params_obpf.attenuation_in_band = 0;
% params_obpf.attenuation_out_band = -20;
tf = opt_tf_obpf(params_obpf,frequency_array);
% Optical bandpass filter transfer function

sig_obpf = opt_filter(sig,tf);
% Filter optical signal to check waveform (eye diagram) 
% Otherwise the signal may be drowned in noise

params_eye.colour_grade = 0;
params_eye.name = ['Optical eye diagram after noise addition'];
meas_eye(sig_obpf,params_eye);
% Plot eye diagram

params_osa.pol = 'both';%'x';%'y';
params_osa.display_interval = [frequency_array(1) frequency_array(end)];
params_osa.resolution_bandwidth = 12.5e9;
params_osa.sensitivity = -100;
params_osa.display = 1;
params_osa.save.txt = 0;
params_osa.save.emf = 0;
params_osa.save.jpg = 0;
params_osa.name = 'Optical spectrum after noise addition';
meas_osa(sig,params_osa); 
text(-100,params_osa.sensitivity + 5,['RB = ' num2str(params_osa.resolution_bandwidth/1.0e9) ' GHz'])
% Check optical spectrum after noise addition


%%
% Extract OSNR

params_osnr.sigfreq = params_tx.emission_frequency - reference_frequency;
params_osnr.meas_bw = 100e9;
params_osnr.meas_offset = 200e9;
params_osnr.noise_bw = 12.5e9;
osnr_retrieved_db = meas_osnr(sig,params_osnr);

fprintf('\n\n%s%i','Number of samples in noise measurement BW: ',floor(params_osnr.meas_bw/df))


fprintf('\n\n%s\t%3.2f%s\n','Target OSNR: ',osnr_target_db,' dB')
fprintf('%s%3.2f%s\n\n','Retreived OSNR: ',osnr_retrieved_db,' dB')

vline((params_osnr.sigfreq + params_osnr.meas_offset + params_osnr.meas_bw/2)/1.0e9,'r--')
vline((params_osnr.sigfreq + params_osnr.meas_offset - params_osnr.meas_bw/2)/1.0e9,'r--')
vline((params_osnr.sigfreq - params_osnr.meas_offset + params_osnr.meas_bw/2)/1.0e9,'r--')
vline((params_osnr.sigfreq - params_osnr.meas_offset - params_osnr.meas_bw/2)/1.0e9,'r--')
vline((params_osnr.sigfreq + params_osnr.meas_bw/2)/1.0e9,'b:')
vline((params_osnr.sigfreq - params_osnr.meas_bw/2)/1.0e9,'b:')
% Plot limits of measurement bandwidths on spectrum 






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
