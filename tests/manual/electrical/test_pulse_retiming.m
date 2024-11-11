% -------------------------------------------------------------------------
% We try different ways of retiming a pulse centered at, say t0, to
% the origin t = 0 (with appropriate periodic boundary conditions).
% This can be done:
% 1. By identifying the pulse maximum, when the pulse is reasonably well
% "peaked". This is all right for standard pulse shapes such as Gaussian
% (1st order), sech, sinc, raised-cosine, root raised-cosine, etc. However,
% this fails for "flat-top" pulses such as rectangular pulses.
% 2. By detrmining the pulse or pulse power (pulse^2) first moment. This
% provides the pulse "center of gravity", which is then retimed to the
% origin.
% 3. By considering the cross-correlation of the pulse with a reference 
% pulse with known timing. 
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
nsamples_per_symbol = 32;
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
% Select pulse type
% -------------------------------------------------------------------------
% pulse_type = 'rrc';
pulse_type = 'rect';
roll_off = 0.1;

pulse_fwhm = 1/symbol_rate/4;
pulse_order = 1;

pulse_position = time_array(nsamples/2);
% Centre position of the pulse.

pulse_position_sample_index = pulse_position/dt + 1;
% Centre position of the pulse, expressed in terms of number of samples.


% -------------------------------------------------------------------------
% Generate isolated pulse
% -------------------------------------------------------------------------
switch pulse_type 
    % Switch over pulse type.
    
    case 'rrc'
        
        sig_pulse = elec_pulse_rrc(time_array,pulse_position,1/symbol_rate,roll_off);         
        
    case 'rc'
        
        sig_pulse = elec_pulse_rc(time_array,pulse_position,1/symbol_rate,roll_off); 
        
    case 'sinc'
        
        sig_pulse = elec_pulse_sinc(time_array,pulse_position,1/symbol_rate);
        
    case 'sech'
     
        sig_pulse = elec_pulse_sech(time_array,pulse_position,pulse_fwhm);
        
    case 'gauss'
        
        sig_pulse = elec_pulse_gauss(time_array,pulse_position,pulse_fwhm,pulse_order);   
        
        
    case 'rect'
        
        nsamples_rect_pulse = 32;
        
        sig_pulse = zeros(1,nsamples);
        
        sig_pulse(pulse_position_sample_index - nsamples_rect_pulse:pulse_position_sample_index + nsamples_rect_pulse) = 1;   
        
        
        
end
% End of switch over pulse type.


figure('Name',[pulse_type ' pulse']);
plot(sig_pulse,'b-')
xlabel('sample')
ylabel('amplitude a.u.')


fprintf('\n\n%s','original pulse:')
fprintf('\n%s\t%i','position (index):',pulse_position_sample_index)
fprintf('\n%s\t%3.15e%s\n','position (time):',pulse_position,' s')



%%
% -------------------------------------------------------------------------
% Recentring the pulse to the origin
% -------------------------------------------------------------------------

% We now attempt to recenter the pulse to another position, for instance to
% the origin, sample index = 1 or pulse position = 0 s.

% -------------------
% Pulse maximum value
% -------------------
[~,sample_index_max] = max(sig_pulse);
% Extract location of the isolated pulse, in term of sample index.
% Will work well in case of non flat-top pulses.

sig_retimed_max = circshift(sig_pulse,-sample_index_max + 1,2);
% Retime the impulse response so that it is centered at zero.

fprintf('\n\n%s','retrieval of pulse position using pulse max:')
fprintf('\n%s\t%i','position (index):',sample_index_max)
fprintf('\n%s\t%3.15e%s\n','position (time):',time_array(sample_index_max),' s')




% -------------------
% Pulse first moment
% -------------------
pulse_moment_time_analog = num_int1d_simpson(time_array.*sig_pulse,dt)./num_int1d_simpson(sig_pulse,dt);
% Pulse first moment (analog).


pulse_moment_sample_digital = round(sum([1:nsamples].*sig_pulse)/sum(sig_pulse));
% Pulse first moment (digital).
% pulse_moment_sample_digital = round(sum([1:nsamples].*sig_pulse.^2)/sum(sig_pulse.^2));
% Pulse first moment (digital) based on pulse power.
% Choose one of the above...

sig_retimed_moment = circshift(sig_pulse,-pulse_moment_sample_digital + 1,2);
% Retime the impulse response so that it is centered at zero.


fprintf('\n\n%s','retrieval of pulse position using pulse 1st moment - digital:')
fprintf('\n%s\t%i','position (index):',pulse_moment_sample_digital)
fprintf('\n%s\t%3.15e%s\n','position (time):',time_array(pulse_moment_sample_digital),' s')

fprintf('\n\n%s','retrieval of pulse position using pulse 1st moment - analog:')
fprintf('\n%s\t%i','position (index):',round(pulse_moment_time_analog/dt + 1))
fprintf('\n%s\t%3.15e%s\n','position (time):',pulse_moment_time_analog,' s')



figure('Name','retimed pulses') 
plot(sig_pulse,'k')
hold on
plot(sig_retimed_max,'b')
plot(sig_retimed_moment,'r--')
xlabel('sample')
ylabel('amplitude (a.u.)')
legend('original','retimed: max pulse value','retimed: 1st moment')



%%
% ---------------------
% Pulse autocorrelation
% ---------------------

% Trying something with autocorrelation...


% We will be using the (adapted) xcorr function from Octave, since we do
% not have access to the Matlab signal processing toolbox.
% Wait a minute. It happens we have a local copy of the Matlab xcorr
% function.
% Let's compare...

fprintf('\n\n\n%s\n','pulse autocorrelation calculation')
fprintf('%s\n','=================================')


fprintf('\n%s\n','octave xcorr')

[ac,lag_ac] = xcorr(sig_pulse);
% Octave xcorr function (made Matlab-compatible)


figure('Name','check pulse autocorrelation')
plot(lag_ac,ac,'b-')
legend('octave xcorr')
xlabel('lag')
ylabel('autocorrelation (a.u.)')





pulse_offset = 33

sig_ref = circshift(sig_pulse,pulse_offset,2);

sig_ref_centre = round(sum([1:nsamples].*sig_ref.^2)/sum(sig_ref.^2));
% Pulse first moment (digital).


figure('Name','reference pulse for autocorrelation')
plot(sig_ref,'r')
xlabel('sample')
ylabel('amplitude (a.u.)')


[ac_pulse,lag_ac_pulse] = xcorr(sig_pulse,sig_ref);
% Cross-correlation between the pulse and the reference pulse.

figure('Name','check pulse autocorrelation')
plot(lag_ac_pulse,ac_pulse,'b-')
xlabel('lag')
ylabel('autocorrelation (a.u.)')



[~,iac_max] = max(ac_pulse);

retrieved_offset = lag_ac_pulse(iac_max)

% We are able to retrieve the delay between a pulse and a delayed version 
% of itself through the autocorrelation.
% Solution with pulse first moment also seems pretty good...






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
