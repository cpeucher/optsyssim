% -------------------------------------------------------------------------
% Test of the char_pulse_rms.m function
%
% We test the calculation of the rms duration of optical pulses in 
% different scenarios:
%  1. Calculation of the rms duration of a first order Gaussian pulse as a
%  function of its FWHM duration.
%  2. Calculation of the second-order dispersion induced broadening of 
%  super-Gaussian pulses.
% 3. Calculation of second- and third-order dispersion induced broadening 
%  of 1st-order Gaussian pulses.
% 
% We take the opportunity to double-check a couple of other functions.
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



        
%%        
% -------------------------------------------------------------------------
% rms duration of a first order Gaussian pulse
% -------------------------------------------------------------------------
gaussian_pulse_fwhm = [10:10:200]*1.0e-12;
% Gaussian pulse FWHM duration, in s

peak_power = 1.0e-3;
chirp = 0;
order = 1;
tc = time_array(length(time_array)/2);
% Gaussian pulse parameters

rms = zeros(1,length(gaussian_pulse_fwhm));

for ipulse=1:length(gaussian_pulse_fwhm)
    
    pulse = opt_pulse_gauss(time_array,order,peak_power,tc,gaussian_pulse_fwhm(ipulse),chirp);
    % Generate single optical pulse
    
    rms(ipulse) = char_pulse_rms(abs(pulse).^2);
    % Calculate rms pulse width numerically

end

fig_name = 'RMS vs FWHM for Gaussian pulse';
hfig = figure('Name',fig_name);
plot(gaussian_pulse_fwhm/1.0e-12,rms/1.0e-12,'Color','b','LineStyle','none','Marker','s','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
hold on
plot(gaussian_pulse_fwhm/1.0e-12,gaussian_pulse_fwhm/1.0e-12/2/sqrt(2*log(2)),'Color','r','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
xlabel('FWHM pulse width (ps)','Interpreter',fig.interpreter)
ylabel('rms pulse width (ps)','Interpreter',fig.interpreter)
legend('numerical','theory','Location','NorthWest','Box','off','Interpreter',fig.interpreter,'FontSize',18)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 18;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
% xlim([0 18])
% ylim([-0.5 0.05])



%%
% -------------------------------------------------------------------------
% Second-order dispersion induced broadening of super-Gaussian pulses
% -------------------------------------------------------------------------

distance_norm = [0:0.2:4];
% Normalised propagation distances (with respect to dispersion length LD) 
% at which the broadening will be calculated

beta2 = -2.17e-26;
% Second order dispersion, in s^2/m

params_pulse.peak_power = 1.0e-3;
params_pulse.chirp = 0;
params_pulse.order = 1;
params_pulse.tc = time_array(length(time_array)/2);
params_pulse.fwhm = 10e-12;
% Gaussian pulse parameters
    
pulse = opt_pulse_gauss(time_array,params_pulse.order,params_pulse.peak_power,...
    params_pulse.tc,params_pulse.fwhm,params_pulse.chirp);
% Generate single optical pulse

sig = struct;
sig.x = pulse;
sig.y = zeros(1,nsamples);
% Convert to optical signal structure

params_scope.visualisers = {'power','chirp'};
params_scope.display_interval = [0 time_array(end)];
params_scope.save.emf = 0;
params_scope.name = 'Check pulse';
meas_scope(sig.x,params_scope); 
% Scope

rms0 = char_pulse_rms(abs(sig.x).^2);
% Calculate rms pulse width at fibre input

pulse_width_t0 = calc_fwhm_t0('gaussian',params_pulse.order,'from_fwhm',params_pulse.fwhm);
% Conversion of pulse duration from FWHM to T0 (half-width at 1/e)

ld = pulse_width_t0^2/abs(beta2);
% Calculation of dispersion length, in m

distance = distance_norm*ld;
% Real distance, in m

params_fibre.nonlinear_coefficient = 0;% nonlinear coefficient, in 1/W/m
params_fibre.dispersion_spec_frequency = reference_frequency;
params_fibre.loss = 0;% loss, in dB/km
params_fibre.beta_coefficients = [0 beta2];
params_fibre.loss_alpha = conv_loss_lin_log(params_fibre.loss);

numparams_fibre.max_step_size = 1;% maximum step size, in m
numparams_fibre.max_phase_shift = 1e-3;% maximum nonlinear phase shift, in radians
% Not relevant here since propagation is linear...

dispersion = conv_disp_d_beta([0 beta2],'from_beta','si','eng',params_fibre.dispersion_spec_frequency);
% Convert dispersion to D, C, S parameters, as required by the
% opt_dispersion function

params_dispersion.dispersion = dispersion(1);
params_dispersion.dispersion_slope = dispersion(2);
params_dispersion.dispersion_curvature = dispersion(3);
params_dispersion.dispersion_spec_frequency = reference_frequency;


rms_nlse = zeros(1,length(distance));
rms_disp = zeros(1,length(distance));


for idist = 1:length(distance)

    params_fibre.length = distance(idist);
    % Fibre length, in m    

    sig_nlse = opt_nlse_scalar_basic(sig,params_fibre,numparams_fibre);
    % Propagation using NLSE
    
    sig_disp = opt_dispersion(sig,params_dispersion,distance(idist)); 
    % Propagation using passive dispersive element

    rms_nlse(idist) = char_pulse_rms(abs(sig_nlse.x).^2);
    rms_disp(idist) = char_pulse_rms(abs(sig_disp.x).^2);
    % Calculate rms pulse duration
    
end
% End of loop over fibre length


bbroadening = sqrt(1 + distance_norm.^2);
% Theoretical broadening for an unchirped Gaussian pulse

rms_broadening_nlse = rms_nlse/rms0;
rms_broadening_disp = rms_disp/rms0;
% Numerically-obtained rms broadening

[rms_broadening_analytical,ld_func] = calc_disp_broadening_gauss_m(distance,beta2,params_pulse);
% We also check the calc_disp_broadening_gauss_m function...

fig_name = 'Broadening of 1st order Gaussian pulse under second order dispersion';
hfig = figure('Name',fig_name);
h1 = plot(distance_norm,rms_broadening_analytical,'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on
h2 = plot(distance_norm,bbroadening,'Color','r','LineStyle','--','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
h3 = plot(distance_norm,rms_broadening_nlse,'Color','b','LineStyle','none','Marker','o','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
h4 = plot(distance_norm,rms_broadening_disp,'Color','k','LineStyle','none','Marker','s','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');

xlabel('normalised distance $z/L_D$','Interpreter',fig.interpreter)
ylabel('rms broadening factor','Interpreter',fig.interpreter)
legend([h1 h2 h3 h4],{'analytical: general expression','analytical: 1st order','numerical: NLSE','numerical: dispersive element'},'Location','NorthWest','Box','off','Interpreter',fig.interpreter,'FontSize',18)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 18;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
% xlim([0 18])
% ylim([-0.5 0.05]) 





%%
% -------------------------------------------------------------------------
% Compare the 2 analytical expressions:
% 1. Broadening of a mth order gaussian pulse under 2nd order dispersion
% (here m=1)
% 2. Broadening of a 1st order chirped Gaussian pulse under 2nd and 3rd 
% order dispersion (here \beta_3 = 0 and C = 0).
% -------------------------------------------------------------------------
params_pulse.peak_power = 1.0e-3;
params_pulse.chirp = 0;
params_pulse.order = 1;
params_pulse.tc = time_array(length(time_array)/2);
params_pulse.fwhm = 10e-12;
% Gaussian pulse parameters

pulse_width_t0 = calc_fwhm_t0('gaussian',params_pulse.order,'from_fwhm',params_pulse.fwhm);
% Conversion of pulse duration from FWHM to T0 (half-width at 1/e)

beta2 = -3.41e-27; 
% beta_2 in s^2/m

ld = pulse_width_t0^2/abs(beta2);
% Dispersion lengths linked to 2nd order dispersion, in m

distance = [0:0.1:5]*ld;

[rms_broadening_analytical_m,ld_func_m] = calc_disp_broadening_gauss_m(distance,beta2,params_pulse);
[rms_broadening_analytical_1,ld_func_1,ldp_func_1] = calc_disp_broadening_gauss_1(distance,beta2,0,params_pulse);
% Pulse rms broadening factor
% We set beta3 to 0 in the call to calc_disp_broadening_gauss_1

fig_name = 'Compare 2 analytical expressions for beta_3 = 0 and m = 1';
hfig = figure('Name',fig_name);
h1 = plot(distance/ld,rms_broadening_analytical_m,'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on
h2 = plot(distance/ld,rms_broadening_analytical_1,'Color','r','LineStyle','--','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');

xlabel('normalised distance $z/L_D$','Interpreter',fig.interpreter)
ylabel('rms broadening factor','Interpreter',fig.interpreter)
legend([h1 h2],{'analytical: order $m$, $\beta_2$','analytical: order 1, $\beta_2$ and $\beta_3$'},'Location','NorthWest','Box','off','Interpreter',fig.interpreter,'FontSize',18)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 18;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
% xlim([0 18])
% ylim([-0.5 0.05])



 

%%
% -------------------------------------------------------------------------
% Comparison between analytical and numerical estimations of broadening for
% a 1st order Gaussian pulse under 2nd and 3rd order dispersion
% -------------------------------------------------------------------------
params_pulse.peak_power = 1.0e-3;
params_pulse.chirp = 0;
params_pulse.order = 1;
params_pulse.tc = time_array(length(time_array)/2);
params_pulse.fwhm = 10e-12;
% Gaussian pulse parameters

pulse_width_t0 = calc_fwhm_t0('gaussian',params_pulse.order,'from_fwhm',params_pulse.fwhm);
% Conversion of pulse duration from FWHM to T0 (half-width at 1/e)

beta2 = -3.41e-27; 
% beta_2 in s^2/m
beta3 = 4.09e-38; 
% beta_3 in s^3/m

ld = pulse_width_t0^2/abs(beta2);
ldp = pulse_width_t0^3/abs(beta3);
% Calculation of dispersion lengths linked to 2nd and 3rd order dispersion,
% in m

distance_norm = [0:0.2:4];
% Distance normalised to the third order dispersion length

pulse = opt_pulse_gauss(time_array,params_pulse.order,params_pulse.peak_power,...
    params_pulse.tc,params_pulse.fwhm,params_pulse.chirp);
% Generate single optical pulse

sig = struct;
sig.x = pulse;
sig.y = zeros(1,nsamples);
% Convert to optical signal structure


params_fibre.nonlinear_coefficient = 0;% nonlinear coefficient, in 1/W/m
params_fibre.dispersion_spec_frequency = reference_frequency;
params_fibre.loss = 0;% loss, in dB/km
params_fibre.beta_coefficients = [0 beta2 beta3];
params_fibre.loss_alpha = conv_loss_lin_log(params_fibre.loss);
% Fibre parameters

numparams_fibre.max_step_size = 1;% maximum step size, in m
numparams_fibre.max_phase_shift = 1e-3;% maximum nonlinear phase shift, in radians
% Not relevant here since propagation is linear...

dispersion = conv_disp_d_beta([0 beta2 beta3],'from_beta','si','eng',params_fibre.dispersion_spec_frequency);
% Convert dispersion to D, C, S parameters, as required by the
% opt_dispersion function

params_dispersion.dispersion = dispersion(1);
params_dispersion.dispersion_slope = dispersion(2);
params_dispersion.dispersion_curvature = dispersion(3);
params_dispersion.dispersion_spec_frequency = reference_frequency;

rms0 = char_pulse_rms(abs(sig.x).^2);
% Calculate rms pulse width at fibre input

distance = distance_norm*ldp;
% Real distance, in m.

rms_nlse = zeros(1,length(distance));
rms_disp = zeros(1,length(distance));


for idist=1:length(distance)
    
    params_fibre.length = distance(idist);
    % Fibre length, in m    

    sig_nlse = opt_nlse_scalar_basic(sig,params_fibre,numparams_fibre);
    % Propagation using NLSE
    
    sig_disp = opt_dispersion(sig,params_dispersion,distance(idist)); 
    % Propagation using passive dispersive element

    rms_nlse(idist) = char_pulse_rms(abs(sig_nlse.x).^2);
    rms_disp(idist) = char_pulse_rms(abs(sig_disp.x).^2);
    % Calculate rms pulse duration
    
    
end
% End of loop over fibre length.


rms_broadening_nlse = rms_nlse/rms0;
rms_broadening_disp = rms_disp/rms0;
% Numerically-obtained rms broadening

[rms_broadening_analytical,ld_func,ldp_func] = calc_disp_broadening_gauss_1(distance,beta2,beta3,params_pulse);
% We also check the calc_disp_broadening_gauss_m function...

fig_name = 'Broadening of 1st order chirped Gaussian pulse under second and third order dispersion';
hfig = figure('Name',fig_name);
h1 = plot(distance_norm,rms_broadening_analytical,'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on
h2 = plot(distance_norm,rms_broadening_nlse,'Color','b','LineStyle','none','Marker','o','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
h3 = plot(distance_norm,rms_broadening_disp,'Color','k','LineStyle','none','Marker','s','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');

xlabel('normalised distance $z/L_D$','Interpreter',fig.interpreter)
ylabel('rms broadening factor','Interpreter',fig.interpreter)
legend([h1 h2 h3],{'analytical','numerical: NLSE','numerical: dispersive element'},'Location','NorthWest','Box','off','Interpreter',fig.interpreter,'FontSize',18)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 18;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
% xlim([0 18])
% ylim([-0.5 0.05]) 





% We also check the broadened pulse after z = 4*ldp

fig_name = 'Waveform of 1st order chirped Gaussian pulse under second and third order dispersion';
hfig = figure('Name',fig_name);
h1 = plot((time_array - time_array(nsamples/2))/1.0e-12,abs(sig.x).^2,'Color','k','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on
h2 = plot((time_array - time_array(nsamples/2))/1.0e-12,abs(sig_nlse.x).^2,'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
h3 = plot((time_array - time_array(nsamples/2))/1.0e-12,abs(sig_disp.x).^2,'Color','r','LineStyle','--','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
xlabel('time (ps)','Interpreter',fig.interpreter)
ylabel('power (W) ','Interpreter',fig.interpreter)
legend([h1 h2 h3],{'input pulse','NLSE','dispersive element'},'Location','NorthEast','Box','off','Interpreter',fig.interpreter,'FontSize',18)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 18;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
xlim([-50 150])
% ylim([-0.5 0.05])


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
