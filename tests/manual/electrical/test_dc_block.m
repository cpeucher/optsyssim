% -------------------------------------------------------------------------
% Test of elec_dc_block function
%
% We test the estimation of the DC level depending on the number of bits
% (when some imbalance exists between the number of marks and spaces).
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
% reference_frequency = 193.1e12;
% nsamples_per_symbol = 128;
% nsymbols = 128;
% symbol_rate = 40e9;
% 
% nsamples = nsamples_per_symbol*nsymbols;
% sample_rate = nsamples_per_symbol*symbol_rate;
% [time_array,dt,frequency_array,df] = core_create_time_axis(nsamples_per_symbol,nsymbols,symbol_rate);        
        
  


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
% 
% -------------------------------------------------------------------------

reference_frequency = 193.1e12;
nsamples_per_symbol = 128;
symbol_rate = 40e9;


% -------------------------------------------------------------------------
% 128 symbols, positive DC value
% -------------------------------------------------------------------------

nsymbols = 128;

nsamples = nsamples_per_symbol*nsymbols;
sample_rate = nsamples_per_symbol*symbol_rate;
[time_array,dt,frequency_array,df] = core_create_time_axis(nsamples_per_symbol,nsymbols,symbol_rate);   

params_prbs.type = 'shift_register';%'de_bruijn';
params_prbs.order = 7;
params_prbs.poly = [7 6 0];%[7 3 0];[7 1 0];
params_prbs.seed = [1 1 1 0 1 1 0];

rise_time = 1/symbol_rate/4;  
% NRZ signal rise time, in s.

bit_pattern = logical_prbs(params_prbs);
bit_pattern = logical_adapt_binary_sequence(bit_pattern,nsymbols);
% Generate random sequence and adapt length to simulation time window

sig = elec_pulse_sequence_nrz(bit_pattern,rise_time);
% Encode the binary bit pattern onto an NRZ signal


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
params_eye.name = 'Original signal';
meas_eye(sig,params_eye);
% Plot eye diagram

[sig_dcb, dc_dcb] = elec_dc_block(sig);

dc_mean = mean(sig);



params_eye.name = 'After DC block';
meas_eye(sig_dcb,params_eye);
% Plot eye diagram


fprintf('\n%s\t%i\n','DC level for Nsymbols =',nsymbols)
fprintf('%s\n','==========================')
fprintf('%s\t%.6f\n','DC level (freq):',dc_dcb)
fprintf('%s\t%6f\n\n\n','DC level (time):',dc_mean)


% -------------------------------------------------------------------------
% 16 symbols, positive DC value
% -------------------------------------------------------------------------
nsymbols = 16;


nsamples = nsamples_per_symbol*nsymbols;
sample_rate = nsamples_per_symbol*symbol_rate;
[time_array,dt,frequency_array,df] = core_create_time_axis(nsamples_per_symbol,nsymbols,symbol_rate);   

params_prbs.type = 'shift_register';%'de_bruijn';
params_prbs.order = 7;
params_prbs.poly = [7 6 0];%[7 3 0];[7 1 0];
params_prbs.seed = [1 1 1 0 1 1 0];

rise_time = 1/symbol_rate/4;  
% NRZ signal rise time, in s.

bit_pattern = logical_prbs(params_prbs);
bit_pattern = logical_adapt_binary_sequence(bit_pattern,nsymbols);
% Generate random sequence and adapt length to simulation time window

sig = elec_pulse_sequence_nrz(bit_pattern,rise_time);
% Encode the binary bit pattern onto an NRZ signal

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
params_eye.name = 'Original signal';
meas_eye(sig,params_eye);
% Plot eye diagram

[sig_dcb, dc_dcb] = elec_dc_block(sig);

dc_mean = mean(sig);



params_eye.name = 'After DC block';
meas_eye(sig_dcb,params_eye);
% Plot eye diagram


fprintf('\n%s\t%i\n','DC level for Nsymbols =',nsymbols)
fprintf('%s\n','==========================')
fprintf('%s\t%.6f\n','DC level (freq):',dc_dcb)
fprintf('%s\t%6f\n','DC level (time):',dc_mean)


% -------------------------------------------------------------------------
% 128 symbols, negative DC value
% -------------------------------------------------------------------------

nsymbols = 128;

nsamples = nsamples_per_symbol*nsymbols;
sample_rate = nsamples_per_symbol*symbol_rate;
[time_array,dt,frequency_array,df] = core_create_time_axis(nsamples_per_symbol,nsymbols,symbol_rate);   

params_prbs.type = 'shift_register';%'de_bruijn';
params_prbs.order = 7;
params_prbs.poly = [7 6 0];%[7 3 0];[7 1 0];
params_prbs.seed = [1 1 1 0 1 1 0];

rise_time = 1/symbol_rate/4;  
% NRZ signal rise time, in s.

bit_pattern = logical_prbs(params_prbs);
bit_pattern = logical_adapt_binary_sequence(bit_pattern,nsymbols);
% Generate random sequence and adapt length to simulation time window

sig = -elec_pulse_sequence_nrz(bit_pattern,rise_time) - 3;
% Encode the binary bit pattern onto an NRZ signal


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
params_eye.name = 'Original signal';
meas_eye(sig,params_eye);
% Plot eye diagram

[sig_dcb, dc_dcb] = elec_dc_block(sig);

dc_mean = mean(sig);



params_eye.name = 'After DC block';
meas_eye(sig_dcb,params_eye);
% Plot eye diagram


fprintf('\n%s\t%i\n','DC level for Nsymbols =',nsymbols)
fprintf('%s\n','==========================')
fprintf('%s\t%.6f\n','DC level (freq):',dc_dcb)
fprintf('%s\t%6f\n\n\n','DC level (time):',dc_mean)






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
