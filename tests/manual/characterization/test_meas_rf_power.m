% -------------------------------------------------------------------------
% Test of meas_rf_power function
%
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
nsymbols = 128*16;
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
% Generate sinusoidal signal
% -------------------------------------------------------------------------
params_rf.frequency = symbol_rate/4;
params_rf.phase = 0;
params_rf.vpp = 2;
params_rf.vdc = 1;
sig = elec_sinusoidal(params_rf);
% Generate electrical sinusoidal signal

load_resistance = 38;
% Load resistance of the analyzers.

rf_power_sine = (params_rf.vpp/2)^2/2/load_resistance;
% Signal RF power, in the load resistance, in W.

dc_power_sine = params_rf.vdc^2/load_resistance;
% Signal DC power, in the load resistance, in W.

fprintf('\n\n%s%3.3f%s\n','Signal peak-to-peak voltage: ',params_rf.vpp,' V');
fprintf('%s%3.3f%s\n','Load resistance: ',load_resistance,' ohm');
fprintf('%s%3.3f%s\n','Signal RF power: ',rf_power_sine,' W');
fprintf('%s%3.3f%s\n','Signal RF power: ',10*log10(rf_power_sine/1.0e-3),' dBm');
fprintf('%s%3.3f%s\n','Signal DC power: ',dc_power_sine,' W');
fprintf('%s%3.3f%s\n','Signal DC power: ',10*log10(dc_power_sine/1.0e-3),' dBm');

xf = fftshift(fft(sig)/nsamples);
% Frequency spectrum of the electrical signal
xf = xf/sqrt(load_resistance);
% Taking the load resistance into account

figure('Name','Double-sided spectrum')
plot(frequency_array,abs(xf).^2,'b')
xlim([-2 2]*params_rf.frequency)
xlabel('frequency (Hz)')
ylabel('power (W)')

fprintf('\n\n%s%3.6f%s\n','Expected power per line in double-sided spectrum: ',rf_power_sine/2,' W');
fprintf('%s%3.3f%s\n','Expected power per line in double-sided spectrum: ',10*log10(rf_power_sine/2/1.0e-3),' dBm');
fprintf('%s%3.3f%s\n','Expected power per line in single-sided spectrum: ',10*log10(rf_power_sine/1.0e-3),' dBm');


params_scope.visualisers = {'amplitude'};
params_scope.display_interval = [0 5/params_rf.frequency];
params_scope.save.emf = 0;
params_scope.name = 'Signal Waveform';
meas_scope(sig,params_scope); 
% Scope

params_esa.display_interval = [0 params_rf.frequency*10];
params_esa.resolution_bandwidth = 0;
params_esa.input_impedance = load_resistance;
params_esa.display = 1;
params_esa.save.txt = 0;
params_esa.save.emf = 0;
params_esa.name = ['RF spectrum - RBW= ' num2str(params_esa.resolution_bandwidth/1.0e6) ' MHz'];
spectrum = meas_esa(sig,params_esa);
% Display RF spectrum in 0 Hz resolution bandwidth

params_esa.display_interval = [0 params_rf.frequency*10];
params_esa.resolution_bandwidth = 9e9;
params_esa.input_impedance = load_resistance;
params_esa.display = 1;
params_esa.save.txt = 0;
params_esa.save.emf = 0;
params_esa.name = ['RF spectrum - RBW= ' num2str(params_esa.resolution_bandwidth/1.0e6) ' MHz'];
spectrum = meas_esa(sig,params_esa);
% Display RF spectrum in 9 GHz resolution bandwidth


params_rf_pwm.centre_frequency = params_rf.frequency;
params_rf_pwm.bandwidth = 4*df;
params_rf_pwm.input_impedance = load_resistance;
[rf_power,psd] = meas_rf_power(sig,params_rf_pwm);
% Measure RF power

fprintf('\n%s%3.3f%s\n','RF power returned from meas_rf_power module: ',10*log10(rf_power/1.0e-3),' dBm');
fprintf('%s%3.6f%s\n','RF power returned from meas_rf_power module: ',rf_power,' W');


params_rf_pwm.centre_frequency = df;
params_rf_pwm.bandwidth = 2*df;
[rf_power,psd] = meas_rf_power(sig,params_rf_pwm);
% Measure RF power

fprintf('\n%s%3.3f%s\n','DC power returned from meas_rf_power module: ',10*log10(rf_power/1.0e-3),' dBm');
fprintf('%s%3.6f%s\n','DC power returned from meas_rf_power module: ',rf_power,' W');


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
