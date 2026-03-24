% -------------------------------------------------------------------------
% We check that the defintion of the transfer function of an electrical 
% filter, as obtained from elec_tf_elpf is correct with respect to the 
% phase and group delay.
%
% We calculate numerically the phase / group delay and make a quick (and
% not robust...) estimation of the phase delay by filtering a sine wave 
% and observing the delay between the output and input (therefore phase
% delay).
%
% Since 
% \tau_g(w) = \tau_phi(w) + w d\tau_phi(w)/dw
% if the variations of the phase delay are correct, so should it be for the
% group delay.
%
% As stated, the estimation of the phase delay is not robust, but just good
% enough for the sake of the present verification of the transfer function
% returned by elec_tf_elpf.
%
% The calculations are performed here for an RC filter.
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
% Distinguishable colors for graphs
% Available at:
% https://se.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap


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
% nsymbols = 8;
% symbol_rate = 10e9;
% 
% nsamples = nsamples_per_symbol*nsymbols;
% sample_rate = nsamples_per_symbol*symbol_rate;
% [time_array,dt,frequency_array,df] = core_create_time_axis(nsamples_per_symbol,nsymbols,symbol_rate);       
        
  
% reference_frequency = 193.1e12;
% df = 1e6;
% nsamples = 2^16;
% 
% sample_rate = nsamples*df;
% dt = 1/sample_rate;
% time_array = (0:nsamples-1)*dt;
% frequency_array = (-nsamples/2:nsamples/2-1)*df;


reference_frequency = 193.1e12;
dt = 1e-12;
nsamples = 2^15;

sample_rate = 1/dt;
df = sample_rate/nsamples;
time_array = (0:nsamples-1)*dt;
frequency_array = (-nsamples/2:nsamples/2-1)*df;

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
% Filter parameters
% -------------------------------------------------------------------------
params_elpf.type = 'rc';%'bessel';%'butterworth';'gaussian;'none';'rc';'rectangular';'raised_cosine';'root_raised_cosine';
% params_elpf.order = 4;
params_elpf.f3dB = 1e9;
% params_elpf.roll_off = 0.6;            % for 'raised_cosine' and 'root_raised_cosine'
% params_elpf.symbol_rate = symbol_rate; % for 'raised_cosine' and 'root_raised_cosine'



% -------------------------------------------------------------------------
% Transfer function analysis
% -------------------------------------------------------------------------
[freq_tf,magnitude,phase,phase_delay,freq_delay,group_delay] = elec_tf_analysis([0 20e9],params_elpf,10000,1);

fig_name = [file_name_core_figure '_' params_elpf.type '_magnitude'];
hfig = figure('Name',fig_name);
semilogx(freq_tf,20*log10(magnitude),'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
xlabel('frequency (Hz)','Interpreter',fig.interpreter)
ylabel('transmission (dB)','Interpreter',fig.interpreter)
legend('','Location','NorthEast','Box','off','Interpreter',fig.interpreter,'FontSize',fig.font_size)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = fig.font_size;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
xticks([1e6 1e7 1e8 1e9 1e10])
% xlim([0 18])
% ylim([-0.5 0.05])

fig_name = [file_name_core_figure '_' params_elpf.type '_phase_delay'];
hfig = figure('Name',fig_name);
semilogx(freq_tf,phase_delay/1.0e-12,'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
xlabel('frequency (Hz)','Interpreter',fig.interpreter)
ylabel('phase delay (ps)','Interpreter',fig.interpreter)
legend('','Location','NorthEast','Box','off','Interpreter',fig.interpreter,'FontSize',fig.font_size)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = fig.font_size;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
xticks([1e6 1e7 1e8 1e9 1e10])
% xlim([0 18])
% ylim([-0.5 0.05])


fig_name = [file_name_core_figure '_' params_elpf.type '_group_delay'];
hfig = figure('Name',fig_name);
semilogx(freq_delay,group_delay/1.0e-12,'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
xlabel('frequency (Hz)','Interpreter',fig.interpreter)
ylabel('group delay (ps)','Interpreter',fig.interpreter)
legend('','Location','NorthEast','Box','off','Interpreter',fig.interpreter,'FontSize',fig.font_size)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = fig.font_size;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
xticks([1e6 1e7 1e8 1e9 1e10])
% xlim([0 18])
% ylim([-0.5 0.05])



%%
% -------------------------------------------------------------------------
% Check phase delay as a function of frequency
% This is done:
% 1. From the calculated transfer function, by applying 
%    \tau_phi = -phi(w)/w
%    where phi is the argument of the transfer function.
% 2. Numerically, by filtering a sine wave and estimating the delay between
%    the filter output and the stimulus (not robust)
% 3. Theoretically
%    For an RC filter:
%    \tau_phi = arctan(w/w_3dB)/w
% -------------------------------------------------------------------------


nsamples_per_period_min = 32;
% Minimum acceptable number of samples per period of the sinusoidal signal

fmax = 1/nsamples_per_period_min/dt;
% Corresponding maximum frequency



freq_stimulus = [1:2:256]*df;
% Frequencies of the RF signal, in Hz

phase_delay_retrieved = zeros(1,length(freq_stimulus));
% Preinitialisation

for ii = 1:length(freq_stimulus)

params_rf.frequency = freq_stimulus(ii);
params_rf.phase = 0;
params_rf.vpp = 1.0;
params_rf.vdc = 0;
sig_s = elec_sinusoidal(params_rf); 
% Stimulus

sig_f = elec_elpf(sig_s,params_elpf);
% Filtered signal

sig_f = sig_f/(max(sig_f) - min(sig_f));
% Normalise to [-0.5 0.5]

sig_1 = sig_s(time_array<= 1/params_rf.frequency/2);
sig_2 = sig_f(time_array<= 1/params_rf.frequency/2);
time_restricted = time_array(time_array<= 1/params_rf.frequency/2);
% Restrict to one half period

c1 = 0 >= sig_1;
c2 = 0 >= sig_2;
% Find zero-crossing point. Assumes signal is "sufficiently" sampled.

s1 = diff(c1);
s2 = diff(c2);
% We could also determine whether crossing arises on a positive or negative
% to make the algorithm robust. 
% We just want to perform a quick check right now. Not necessary at this
% point.

k1 = find(c1,1);
k2 = find(c2,1);
phase_delay_retrieved(ii) = (k2 -k1)*dt;

if ii == 1

    figure('Name',['Stimulus and filtered signal for f =' num2str(ii*df) ' Hz'])
    plot(time_restricted,sig_1,'b')
    hold on
    plot(time_restricted,sig_2,'r--')
    xlabel('time (s)')
    ylabel('amplitude (a.u.)')
    legend('stimulus','filtered')

end

if ii == length(freq_stimulus)

    figure('Name',['Stimulus and filtered signal for f =' num2str(ii*df) ' Hz'])
    plot(time_restricted,sig_1,'b')
    hold on
    plot(time_restricted,sig_2,'r--')
    xlabel('time (s)')
    ylabel('amplitude (a.u.)')
    legend('stimulus','filtered')

end

end


freq = frequency_array(frequency_array>= 0);
% Positive frequencies
phase_delay_theory = atan(freq/params_elpf.f3dB)./freq/2/pi;
% Theoretical phase delay for an RC filter





fig_name = [file_name_core_figure '_' params_elpf.type '_phase_delay_retrieved'];
hfig = figure('Name',fig_name);
semilogx(freq_tf,phase_delay/1.0e-12,'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
hold on
semilogx(freq_stimulus,phase_delay_retrieved/1.0e-12,'Color','r','LineStyle','none','Marker','o','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
semilogx(freq,phase_delay_theory/1.0e-12,'Color','r','LineStyle','--','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')

xlabel('frequency (Hz)','Interpreter',fig.interpreter)
ylabel('phase delay (ps)','Interpreter',fig.interpreter)
legend('retrieved from elec_tf_elpf','retrieved through simulations','analytical','Location','SouthWest','Box','off','Interpreter',fig.interpreter,'FontSize',fig.font_size)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = fig.font_size;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
xticks([1e6 1e7 1e8 1e9 1e10])
% xlim([0 18])
 ylim([0 160])
grid on



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
