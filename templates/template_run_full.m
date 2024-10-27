% -------------------------------------------------------------------------
% 
%
% 
% 2024-08-01
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
time_stamp = datestr(datetime('now','TimeZone','Z'),'yyyymmddThhMMSSZ');

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

% -------------------------------------------------------------------------
% PRBS
% -------------------------------------------------------------------------
params_prbs.type = 'shift_register';%'de_bruijn';
params_prbs.order = 7;
params_prbs.poly = [7 6 0];%[7 3 0];[7 1 0];
params_prbs.seed = [1 1 1 0 1 1 0];

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
% Generate random sequence and adapt length to simulation time window. 


nrz_data_sig = elec_pulse_sequence_nrz(bit_pattern,rise_time);
% Encode the binary bit pattern onto an NRZ signal.


params_rz.duty_cycle = 0.5;
params_rz.pulse_position = 0.5;
params_rz.rise_time = 1/symbol_rate/8;
params_rz.amplitude_mode = 'normalised';
rz_data_sig = elec_coder_rz(bit_pattern,nsamples_per_symbol,...
          params_rz.duty_cycle,params_rz.pulse_position,...
          params_rz.rise_time,params_rz.amplitude_mode);
      
      
params_scope.type = 'elec';
params_scope.pol = 'x';%'y','both';
params_scope.visualiser_type = 'power';%'power_phase';'power_chirp';'power_phase_chirp';
params_scope.display_interval = [0 time_array(end)];
params_scope.save.ascii = 0;
params_scope.save.emf = 0;
params_scope.name = 'Electrical signal';
meas_scope(rz_data_sig,params_scope);
% Oscilloscope

params_tx.type = 'ook_nrz';
params_tx.emission_frequency = reference_frequency;
params_tx.power = 1.0e-3;
params_tx.linewidth = 0;
params_tx.bit_pattern = bit_pattern;
params_tx.symbol_rate = symbol_rate;
params_tx.rise_time = 1/symbol_rate/4;
sig = tx(params_tx);
% Optical transmitter.




params_osa.pol = 'x';%'y','both';
params_osa.display_interval = [-6*symbol_rate 6*symbol_rate];
params_osa.resolution_bandwidth = 12.5e9;
params_osa.sensitivity = -60;
params_osa.display = 1;
params_osa.save.ascii = 0;
params_osa.save.emf = 0;
params_osa.save.jpg = 0;
params_osa.name = 'Optical spectrum';
meas_osa(sig,params_osa);
% Optical spectrum analyser.


params_scope.type = 'opt';%'elec';
params_scope.pol = 'x';%'y','both';
params_scope.visualiser_type = 'power';%'power_phase';'power_chirp';'power_phase_chirp';
params_scope.display_interval = [0 time_array(end)];
params_scope.save.ascii = 0;
params_scope.save.emf = 0;
params_scope.name = 'Optical signal';
meas_scope(sig,params_scope);
% Oscilloscope



params_pin.pd.responsivity = 1;
params_pin.pd.thermal_noise_density = 10e-12;
params_pin.pd.shot_noise = 0;
params_pin.pd.dark_current = 0;
params_pin.elpf.type = 'bessel';%'none';'gaussian';'rc';'rectangular';
params_pin.elpf.order = 4;
params_pin.elpf.f3dB = 0.7*symbol_rate;
sig = rx_pin(sig,params_pin);
% PIN receiver.
 

params_eye.pol = 'x';%'y','both';
params_eye.eye_number = 2;
params_eye.samples_per_symbol = nsamples_per_symbol;
params_eye.save.ascii = 0;
params_eye.save.emf = 0;
params_eye.save.jpg = 0;
params_eye.save.vertical_scale_type = 'auto';%'fixed';
params_eye.save.vertical_scale_max = 1.0e-3;
params_eye.save.display_baseline = 1;
params_eye.colour_grade = 0;
params_eye.name = 'Eye diagram';
meas_eye(sig,params_eye);
% Eye diagram

params_esa.display_interval = [0 4*symbol_rate];
params_esa.resolution_bandwidth = 10e9;
params_esa.input_impedance = 1;
params_esa.display = 1;
params_esa.save.ascii = 0;
params_esa.save.emf = 0;
params_esa.name = 'RF spectrum';
spectrum = meas_esa(sig,params_esa);
% RF spectrum analyser

params_rf_pwm.centre_frequency = 60e9;
params_rf_pwm.bandwidth = 10e9;
params_rf_pwm.input_impedance = 1;
[rf_power,psd] = meas_rf_power(sig,params_rf_pwm);
% RF power meter.






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
