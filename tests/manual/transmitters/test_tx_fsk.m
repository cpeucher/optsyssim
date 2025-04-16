% -------------------------------------------------------------------------
% Test of tx_fsk function
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
% Generate and check electrical waveforms
% -------------------------------------------------------------------------
params_prbs.type = 'shift_register';
params_prbs.order = 11;
params_prbs.poly = [11 9 0];%[11 8 5 2 0];
params_prbs.seed = [0 1 1 0 1 1 1 0 0 1 0];
bit_pattern = logical_prbs(params_prbs);
bit_pattern = logical_adapt_binary_sequence(bit_pattern,nsymbols);
% Generate random sequence and adapt length to simulation time window

rise_time = 1/symbol_rate/4;  
% NRZ signal rise time, in s

nrz_data_sig = elec_pulse_sequence_nrz(bit_pattern,rise_time);
% Encode the binary bit pattern onto an NRZ signal.


params_scope.visualisers = {'amplitude'};
params_scope.display_interval = [0 32/symbol_rate];
params_scope.save.emf = 0;
params_scope.name = 'Original NRZ data signal';
meas_scope(nrz_data_sig,params_scope); 
% Scope

%%
% -------------------------------------------------------------------------
% Generate CPFSK signal
% -------------------------------------------------------------------------
params_tx.rise_time = 1/symbol_rate/5;
params_tx.emission_frequency = reference_frequency;
params_tx.symbol_rate = symbol_rate;
params_tx.power = 1.0e-3;
params_tx.bt = 0.2;           % For GMSK only.
params_tx.tone_spacing = 50e9;
params_tx.type = 'cpfsk';
sig = tx_fsk(bit_pattern,params_tx);

figure('Name','CPFSK complex envelope')
plot(time_array*symbol_rate,real(sig.x),'b-');
hold on
plot(time_array*symbol_rate,imag(sig.x),'r--');
legend('real part','imaginary part')
xlim([10 15]);
xlabel('bit number');
ylabel('field envelope');

params_scope.visualisers = {'power','chirp'};
params_scope.display_interval = [0 32/symbol_rate];
params_scope.save.emf = 0;
params_scope.name = 'CPFSK signal';
meas_scope(sig.x,params_scope); 
% Scope

params_osa.pol = 'x';%'y';'both';
params_osa.display_interval = [frequency_array(1) frequency_array(end)];
params_osa.resolution_bandwidth = 12.5e9;
params_osa.sensitivity = -80;
params_osa.display = 1;
params_osa.save.txt = 0;
params_osa.save.emf = 0;
params_osa.save.jpg = 0;
params_osa.name = 'CPFSK spectrum';
meas_osa(sig,params_osa); 
% Optical spectrum analyser

%%
% -------------------------------------------------------------------------
% Generate FSK signal
% -------------------------------------------------------------------------
params_tx.type = 'fsk';
params_tx.rise_time = 0;
sig = tx_fsk(bit_pattern,params_tx);

figure('Name','FSK complex envelope')
plot(time_array*symbol_rate,real(sig.x),'b-');
hold on
plot(time_array*symbol_rate,imag(sig.x),'r--');
legend('real part','imaginary part')
xlim([10 15]);
xlabel('bit number');
ylabel('field envelope');

params_scope.visualisers = {'power','chirp'};
params_scope.display_interval = [0 32/symbol_rate];
params_scope.save.emf = 0;
params_scope.name = 'FSK signal';
meas_scope(sig.x,params_scope); 
% Scope

params_osa.pol = 'x';%'y';'both';
params_osa.display_interval = [frequency_array(1) frequency_array(end)];
params_osa.resolution_bandwidth = 12.5e9;
params_osa.sensitivity = -80;
params_osa.display = 1;
params_osa.save.txt = 0;
params_osa.save.emf = 0;
params_osa.save.jpg = 0;
params_osa.name = 'FSK spectrum';
meas_osa(sig,params_osa); 
% Optical spectrum analyser



%%
% -------------------------------------------------------------------------
% Comparison of MSK and GMSK with various BT
% -------------------------------------------------------------------------
params_osa.display = 0;
params_osa.resolution_bandwidth = 1e9;

params_tx.type = 'msk';
sig_msk = tx_fsk(bit_pattern,params_tx);
[raw_msk, averaged_msk] = meas_osa(sig_msk,params_osa);

params_tx.type = 'gmsk';
params_tx.bt = 1;
sig_gmsk_bt10 = tx_fsk(bit_pattern,params_tx);
[raw_gmsk_bt10, averaged_gmsk_bt10] = meas_osa(sig_gmsk_bt10,params_osa);

params_tx.type = 'gmsk';
params_tx.bt = 0.7;
sig_gmsk_bt07 = tx_fsk(bit_pattern,params_tx);
[raw_gmsk_bt07, averaged_gmsk_bt07] = meas_osa(sig_gmsk_bt07,params_osa);

params_tx.type = 'gmsk';
params_tx.bt = 0.5;
sig_gmsk_bt05 = tx_fsk(bit_pattern,params_tx);
[raw_gmsk_bt05, averaged_gmsk_bt05] = meas_osa(sig_gmsk_bt05,params_osa);

params_tx.type = 'gmsk';
params_tx.bt = 0.3;
sig_gmsk_bt03 = tx_fsk(bit_pattern,params_tx);
[raw_gmsk_bt03, averaged_gmsk_bt03] = meas_osa(sig_gmsk_bt03,params_osa);

params_tx.type = 'gmsk';
params_tx.bt = 0.2;
sig_gmsk_bt02 = tx_fsk(bit_pattern,params_tx);
[raw_gmsk_bt02, averaged_gmsk_bt02] = meas_osa(sig_gmsk_bt02,params_osa);


figure('Name','Compare MSK spectra')
plot(frequency_array/params_tx.symbol_rate,10*log10(averaged_msk/max(averaged_msk)),'b-');
hold on;
plot(frequency_array/params_tx.symbol_rate,10*log10(averaged_gmsk_bt10/max(averaged_gmsk_bt10)),'r-');
plot(frequency_array/params_tx.symbol_rate,10*log10(averaged_gmsk_bt07/max(averaged_gmsk_bt07)),'c-');
plot(frequency_array/params_tx.symbol_rate,10*log10(averaged_gmsk_bt05/max(averaged_gmsk_bt05)),'g-');
plot(frequency_array/params_tx.symbol_rate,10*log10(averaged_gmsk_bt03/max(averaged_gmsk_bt03)),'m-');
plot(frequency_array/params_tx.symbol_rate,10*log10(averaged_gmsk_bt02/max(averaged_gmsk_bt02)),'k-');
xlim([0 2.5]);
ylim([-120 5]);
xlabel('normalised frequency');
ylabel('normalised power (dB)');
grid on;
h=legend('MSK','GMSK BT = 1','GMSK BT = 0.7','GMSK BT = 0.5','GMSK BT = 0.3','GMSK BT = 0.2');
set(h,'Interpreter','none','Location','SouthWest','Box','on','Orientation','vertical');


figure('Name','Compare MSK and GMSK waveforms')
plot(time_array*symbol_rate,real(sig_msk.x));
hold on;
plot(time_array*symbol_rate,real(sig_gmsk_bt03.x));
xlim([10 20]);
xlabel('bit number');
ylabel('real part of field');
h=legend('MSK','GMSK BT = 0.3');
set(h,'Interpreter','none','Location','NorthWest','Box','off','Orientation','vertical');   
        





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
