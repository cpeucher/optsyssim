% -------------------------------------------------------------------------
% Test of set_er function
%
% We adjust the extinction ratio of an on-off keying optical signal
% generated using the tx module.
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

set(groot, "defaultFigurePosition", [680 458 560 420])
% Temporary fix so that the figures generated in R2025a have the same
% appearance as those generated in R2024b and earlier.


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

line_color = linspecer(6,'qualitative');
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

params_prbs.type = 'shift_register';%'de_bruijn';
params_prbs.order = 7;
params_prbs.poly = [7 6 0];%[7 3 0];[7 1 0];
params_prbs.seed = [1 1 1 0 1 1 0];
bit_pattern = logical_prbs(params_prbs);
bit_pattern = logical_adapt_binary_sequence(bit_pattern,nsymbols);
% Generate random sequence and adapt length to simulation time window 


sig_generator = 'tx';
% The tx module, generating a signal having an infinite extinction ratio
% will be used.
% However, this creates an issue in case some linewidth is present. Read
% below.
sig_generator = 'mzm';
% By controlling the MZM, we can generate a signal with finite extinction
% ratio that is furthermore compatible with phase noise.




params_tx.type = 'ook_nrz';
params_tx.emission_frequency = reference_frequency;
params_tx.power = 1.0e-3;
params_tx.linewidth = 1e6;
params_tx.bit_pattern = bit_pattern;
params_tx.symbol_rate = symbol_rate;
params_tx.rise_time = 1/symbol_rate/4;


switch sig_generator

    case 'tx'

        sig_tx = tx(params_tx);
        % The signal generated by the tx module has an infinite extinction ratio.
        % Therefore, even though we imprint some phase noise through the linewidth
        % parameter, the phase is undefined whenever a space ('0') is generated.
        % This is not something that can be recovered when introducing a finite
        % extinction ratio to the signal.

    case 'mzm'

        nrz_data_sig = elec_pulse_sequence_nrz(params_tx.bit_pattern,params_tx.rise_time);
        % Generate NRZ data stream

        sig = opt_laser_cw(params_tx);
        % CW laser

        vpi = 1.0;
        driving_signal_1 = 0.9*vpi/2*(nrz_data_sig - 0.5);
        driving_signal_2 = -0.9*vpi/2*(nrz_data_sig - 0.5);
        % We adjust Vpp to ensure the extinction ratio is not infinite.
        bias_1 = 1.5*vpi;
        bias_2 = 0;
        split_in = 0.5;
        split_out = 0.5;
        loss = 0;
        sig_tx = mod_mzm(sig,driving_signal_1,driving_signal_2,bias_1,bias_2,vpi,split_in,split_out,loss);
        % Data modulator

end


Pav_tx_dBm = 10*log10(char_opt_average_power(sig_tx)/1.0e-3)
% Average power at Tx output, in dBm

ER_tx_dB = 10*log10(max(abs(sig_tx.x).^2)/min(abs(sig_tx.x).^2))
% Extinction ratio at Tx output, in dB


params_eye.pol = 'x';%'y','both';
params_eye.neyes = 2;
params_eye.nsamples_per_symbol = nsamples_per_symbol;
params_eye.save.txt = 0;
params_eye.save.emf = 0;
params_eye.save.jpg = 0;
params_eye.save.vertical_scale_type = 'auto';%'fixed';
params_eye.save.vertical_scale_max = 7.0e-3;
params_eye.save.display_baseline = 1;
params_eye.colour_grade = 0;
params_eye.name = 'Original';
meas_eye(sig_tx,params_eye);
% Eye diagram at initial transmitter output

Pav_target_dBm = 0
% Target average power, in dBm
ER_target_dB = 2
% Target extinction ratio, in dB



sig = set_er(sig_tx,ER_target_dB,Pav_target_dBm);
% Set average power and extinction ratio

Pav_dBm = 10*log10(char_opt_average_power(sig)/1.0e-3)
% Retrieved average power, in dBm

ER_dB = 10*log10(max(abs(sig.x).^2)/min(abs(sig.x).^2))
% Retrieved extinction ratio, in dB

params_eye.name = 'After Pav and ER adjustment';
meas_eye(sig,params_eye);
% Eye diagram

figure('Name','Input and output power')
plot(time_array/1.0e-12,abs(sig_tx.x).^2/1.0e-3,'b-')
hold on
plot(time_array/1.0e-12,abs(sig.x).^2/1.0e-3,'r--')
xlabel('time (ps)')
ylabel('power (mW)')
legend('Original','After Pav and ER adjustment')


figure('Name','Input and output phases')
plot(time_array/1.0e-12,-unwrap(angle(sig_tx.x))/pi,'b-')
hold on
plot(time_array/1.0e-12,-unwrap(angle(sig.x))/pi,'r--')
xlabel('time (ps)')
ylabel('phase (\times \pi rad)')
legend('Original','After Pav and ER adjustment')


















 







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
