% -------------------------------------------------------------------------
% Test of elec_coder_rz
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
nsamples_per_symbol = 16;
nsymbols = 4;
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
% Parameters
% -------------------------------------------------------------------------
duty_cycle = 0.6 %0.3;0.5;1;
% RZ pulse duty cycle
rise_time_fraction = 0.25;
% Rise time, expressed in fraction of the symbol slot

% -------------------------------------------------------------------------
% Alternate sequence
% -------------------------------------------------------------------------

isymbol = (1:1:nsymbols);
% Symbol indices

seq_alternate = mod(isymbol,2);
% Alternate 101010... sequence


sig_nrz = elec_coder_nrz(seq_alternate,nsamples_per_symbol);
% NRZ signal


nsamples_per_pulse = round(duty_cycle*nsamples_per_symbol/2)*2;

duty_cycle_actual = nsamples_per_pulse/nsamples_per_symbol

clock = [ones(1,nsamples_per_pulse) zeros(1,nsamples_per_symbol - nsamples_per_pulse)];
clock = repmat(clock,1,nsymbols);
clock = circshift(clock,[0 -(nsamples_per_pulse - nsamples_per_symbol)/2]);

sig_rz = double(and(sig_nrz,clock));

rise_time = rise_time_fraction*nsamples_per_symbol*dt;

sig_rz_rt = elec_rise_time(sig_rz,rise_time);


figure('Name','signal samples')
subplot(4,1,1)
stem((1:nsamples),sig_nrz,'b-')
xlim([1 nsamples])
subplot(4,1,2)
stem((1:nsamples),clock,'r--^')
xlim([1 nsamples])
subplot(4,1,3)
stem((1:nsamples),sig_rz,'r--^')
xlim([1 nsamples])
subplot(4,1,4)
stem((1:nsamples),sig_rz_rt,'r--^')
xlim([1 nsamples])

% -------------------------------------------------------------------------
% Try with functionalised elec_coder_rz
% -------------------------------------------------------------------------
params_rz.duty_cycle = duty_cycle;
params_rz.rise_time = 1/symbol_rate*rise_time_fraction;
params_rz.normalisation = 0;
sig_rz_alt_func = elec_coder_rz(seq_alternate,nsamples_per_symbol,params_rz); 

figure('Name','Compare functionalised')
plot((1:nsamples),sig_rz_alt_func,'b-')
hold on
stem((1:nsamples),sig_rz_rt,'r^')
xlim([1 nsamples])
legend('elec_coder_rz','direct')
xlabel('sample')
ylabel('amplitude (a.u.)')

% -------------------------------------------------------------------------
% Try with PRBS
% -------------------------------------------------------------------------
nsymbols = 16;
nsamples = nsamples_per_symbol*nsymbols;
[time_array,dt,frequency_array,df] = core_create_time_axis(nsamples_per_symbol,nsymbols,symbol_rate);  
% We increase the number of symbols and redefine the time axis accordingly.

params_prbs.type = 'shift_register';%'de_bruijn';
params_prbs.order = 7;
params_prbs.poly = [7 6 0];%[7 3 0];[7 1 0];
params_prbs.seed = [1 1 1 0 1 1 0];
seq_prbs = logical_prbs(params_prbs);
seq_prbs = logical_adapt_binary_sequence(seq_prbs,nsymbols);

sig_nrz_prbs = elec_coder_nrz(seq_prbs,nsamples_per_symbol);

[sig_rz_prbs,params_rz] = elec_coder_rz(seq_prbs,nsamples_per_symbol,params_rz); 

figure('Name','PRBS')
plot((1:nsamples),sig_nrz_prbs,'r')
hold on
plot((1:nsamples),sig_rz_prbs,'b')
xlim([1 nsamples])
xlabel('sample')
ylabel('amplitude (a.u.)')
legend('NRZ','RZ')
ylim([0 1]*1.1)




















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
