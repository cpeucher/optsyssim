% -------------------------------------------------------------------------
% Test of electrical matched filtering
% 
% We synthesise and characterize the transfer function of a matched filter
% for root raised-cosine pulse (possibly other pulse shapes).
% We apply the matched filter to PAM2, PAM4 and QAM16 modulated signals.
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
% Select pulse type
% -------------------------------------------------------------------------
pulse_type = 'rrc';
roll_off = 0.1;

pulse_fwhm = 1/symbol_rate/4;
pulse_order = 1;


% -------------------------------------------------------------------------
% Generate isolated pulse
% -------------------------------------------------------------------------
switch pulse_type 
    % Switch over pulse type.
    
    case 'rrc'
        
        sig_pulse = elec_pulse_rrc(time_array,time_array(nsamples/2),1/symbol_rate,roll_off);         
        
    case 'rc'
        
        sig_pulse = elec_pulse_rc(time_array,time_array(nsamples/2),1/symbol_rate,roll_off); 
        
    case 'sinc'
        
        sig_pulse = elec_pulse_sinc(time_array,time_array(nsamples/2),1/symbol_rate);
        
    case 'sech'
     
        sig_pulse = elec_pulse_sech(time_array,time_array(nsamples/2),pulse_fwhm);
        
    case 'gaussian'
        
        sig_pulse = elec_pulse_gauss(time_array,time_array(nsamples/2),pulse_fwhm,pulse_order);   
        
        
end
% End of switch over pulse type.


figure('Name',['isolated pulse (lin): ' pulse_type]);
plot(sig_pulse,'b-')
xlabel('sample')
ylabel('amplitude a.u.')


figure('Name',['isolated pulse (log): ' pulse_type]);
plot(10*log10(abs(sig_pulse)),'b-')
xlabel('sample')
ylabel('amplitude a.u.')

% We expect a little bit of Gibbs oscillations.
% We could apply a windowing function (we actually did), but this is not
% strictly necessary unless one really wants to avoid tiny oscillations in
% the matched filter transfer function.



ep_pulse = num_int1d_simpson(abs(sig_pulse).^2,dt)
% Calculate pulse energy.

[~,pulse_location] = max(sig_pulse)
% Extract location of the isolated pulse.



%%
% -------------------------------------------------------------------------
% Calculate and characterise matched electrical filter
% -------------------------------------------------------------------------
hmf = fliplr(sig_pulse(:).');
% Impulse response of the matched filter.

K1 = 1/num_int1d_simpson(hmf,dt)

K2 = 1/num_int1d_simpson(hmf.^2,dt)

symbol_rate



K = 1/num_int1d_simpson(hmf,dt);
% Normalisation so that the transfer function of the matched filter at DC
% is equal to 1.
% Note that normally, we do not care about the normalisation, since the
% signal and the noise are scaled with the same value. 
% We just want to check we understand some of the details right.

hmf = K*hmf;
% This is the normalised impulse response.

figure('Name','matched filter impulse response (normalised so the DC attenuation is 1)') 
plot(hmf,'b')


[~,index_max_hmf] = max(hmf);
hmf = circshift(hmf,-index_max_hmf + 1,2);
% Retime the impulse response so that it is centered at zero.

hold on
plot(hmf,'r')
xlabel('sample')
ylabel('amplitude (a.u.)')
legend('original','retimed at origin')


Hmf = num_ft(hmf,dt);
% Transfer function of matched filter, in increasing frequency order.

figure('Name','matched filter transfer function')
plot(frequency_array,10*log10(abs(Hmf).^2),'b')
xlabel('frequency (Hz)')
ylabel('magnitude (dB)')


H0 = num_int1d_simpson(hmf,dt);
% Expected attenuation of the matched filter at DC.

H0db = 10*log10(abs(H0).^2)
% DC attenuation, in dB.

% The eye opening should be, in the absence of residual ISI (here we are
% concerned with a single pulse, therefore this should correpsonds to the
% centre value of the pulse after matched filtering.

eye_opening = K*ep_pulse
% This the value of the sample of the filtered signal at the pulse maximum. 






%%
% -------------------------------------------------------------------------
% Apply matched filter to isolated pulse
% -------------------------------------------------------------------------
% sig = real(ifft(fft(sig).* fftshift(tf)));
% This is the way we normally deal with an electrical filter.
% tf is defined in increasing frequency order, hence the need to fftshift
% it so that it can be applied to fft(sig), which is in fft order.

sig = ifft(fft(sig_pulse).* fftshift(Hmf));
% Apply matched filter.

figure('Name','single pulse after matched filtering')
plot(sig,'b')
xlabel('sample')
ylabel('pulse amplitude (a.u.)')

[max_sig_mf,index_max_sig_mf] = max(sig)

% We check that the max of the matched filter signal is equal to the
% previously calculated "eye opening", i.e. the value in eye_opening, and
% that the filtered pulse is centered at the same time instant as the
% original pulse, i.e. the sample number stored in "pulse_location".



%%
% -------------------------------------------------------------------------
% Generate PAM2 signal with root-raised cosine pulse and apply matched
% filtering.
% -------------------------------------------------------------------------
data = generate_binary(nsymbols,2);
% Generate random binary data.

symbs_pam2 = data;
% For binary modulation.


params_elecmod.pulse_shape = pulse_type;
params_elecmod.roll_off = roll_off;
params_elecmod.symbol_rate = symbol_rate;
params_elecmod.fwhm = pulse_fwhm;
params_elecmod.order = pulse_order;
params_elecmod.normalisation = 'default';
sig_mod = elec_modulator(symbs_pam2,params_elecmod);
% Modulator

eye = calc_eye(sig_mod,nsamples_per_symbol);
% Calculate eye diagram.

figure('Name',['eye diagram (PAM2): ' pulse_type])
plot(eye,'b','LineWidth',1.5);
xlabel('sample')
ylabel('amplitude (a.u.)')
% Plot eye diagram.

constellation_type = 'plain';%'heat';'cluster';
constellation_name = ['constellation (PAM2): ' pulse_type];
plot_constellation(sig_mod(nsamples_per_symbol/2:nsamples_per_symbol:end),constellation_type,constellation_name,[-5:1:5]);
% Constellation



sig_mod_filt = ifft(fft(sig_mod).* fftshift(Hmf));
% Matched filtering.

eye = calc_eye(sig_mod_filt,nsamples_per_symbol);
% Calculate eye diagram.
figure('Name',['eye diagram (PAM2) after matched filter: ' pulse_type])
plot(eye,'b','LineWidth',1.5);
xlabel('sample')
ylabel('amplitude (a.u.)')
% Plot eye diagram.


constellation_name = ['constellation (PAM2) after matched-filtering: ' pulse_type];
plot_constellation(sig_mod_filt(nsamples_per_symbol/2:nsamples_per_symbol:end),constellation_type,constellation_name,[-5:1:5]);
% Constellation


figure('Name','compare matched-filtered and PAM2 original signals')
plot(sig_mod,'b')
hold on
plot(sig_mod_filt,'r')
legend('modulator output','after matched-filtering')
xlabel('sample')
ylabel('amplitude (a.u.)')


% We observe:
% 1. The pulses are properly timed, i.e. no delay is introduced by the
% matched filter.
% 2. The ISI that was introduced by the root raised-cosine pulse shaping is
% suppressed by the matched filter.



%%
% -------------------------------------------------------------------------
% Generate PAM4 signal with root-raised cosine pulse and apply matched
% filtering.
% -------------------------------------------------------------------------

m = 4;

data = generate_binary(nsymbols,m);
% Generate binary data.

[constellation_gray,~,~] = define_constellation('pam4_gray',m);
% Define the constellation for Gray mapping.

[words_dec,words_bin] = conv_bin2dec(data,log2(m));
% Convert binary data into 2-bit words.

symbs_pam4 = constellation_gray(words_dec + 1);
% Symbols generation

sig_mod = elec_modulator(symbs_pam4,params_elecmod);
% Modulator


eye = calc_eye(sig_mod,nsamples_per_symbol);
% Calculate eye diagram.

figure('Name',['eye diagram (PAM2): ' pulse_type])
plot(eye,'b','LineWidth',1.5);
xlabel('sample')
ylabel('amplitude (a.u.)')
% Plot eye diagram.

constellation_type = 'plain';%'heat';'cluster';
constellation_name = ['constellation (PAM4): ' pulse_type];
plot_constellation(sig_mod(nsamples_per_symbol/2:nsamples_per_symbol:end),constellation_type,constellation_name,[-5:1:5]);
% Constellation



sig_mod_filt = ifft(fft(sig_mod).* fftshift(Hmf));
% Matched filtering.

eye = calc_eye(sig_mod_filt,nsamples_per_symbol);
% Calculate eye diagram.
figure('Name',['eye diagram (PAM4) after matched filter: ' pulse_type])
plot(eye,'b','LineWidth',1.5);
xlabel('sample')
ylabel('amplitude (a.u.)')
% Plot eye diagram.


constellation_name = ['constellation (PAM4) after matched-filtering: ' pulse_type];
plot_constellation(sig_mod_filt(nsamples_per_symbol/2:nsamples_per_symbol:end),constellation_type,constellation_name,[-5:1:5]);
% Constellation


figure('Name','compare matched-filtered and original PAM4 signals')
plot(sig_mod,'b')
hold on
plot(sig_mod_filt,'r')
legend('modulator output','after matched-filtering')
xlabel('sample')
ylabel('amplitude (a.u.)')


% We observe:
% 1. The pulses are properly timed, i.e. no delay is introduced by the
% matched filter.
% 2. The ISI that was introduced by the root raised-cosine pulse shaping is
% suppressed by the matched filter.



%%
% -------------------------------------------------------------------------
% Generate QAM16 signal with root-raised cosine pulse and apply matched
% filtering.
% -------------------------------------------------------------------------
m = 16;
% Constellation order. M points with log2(M) even (square constellation)

[constellation,~,~] = define_constellation('qam16_gray',m);
% Define the constellation for Gray mapping.

data = generate_binary(nsymbols,m);
% Generate binary data.

[words_dec,words_bin] = conv_bin2dec(data,log2(m));
% Convert binary data into 2-bit words.

symbs_qam16 = mapping(words_dec,constellation);
% Mapping.

sig_mod = elec_modulator(symbs_qam16,params_elecmod);
% Modulator


eye = calc_eye(sig_mod,nsamples_per_symbol);
% Calculate eye diagram.

figure('Name',['eye diagram (QAM16) - in-phase: ' pulse_type])
plot(real(eye),'b','LineWidth',1.5);
xlabel('sample')
ylabel('amplitude (a.u.)')

figure('Name',['eye diagram (QAM16) - quadrature: ' pulse_type])
plot(imag(eye),'b','LineWidth',1.5);
xlabel('sample')
ylabel('amplitude (a.u.)')

figure('Name',['eye diagram (QAM16) - power: ' pulse_type])
plot(abs(eye).^2,'b','LineWidth',1.5);
xlabel('sample')
ylabel('amplitude (a.u.)')
% Plot eye diagram.

constellation_type = 'plain';%'heat';'cluster';
constellation_name = ['constellation (QAM16): ' pulse_type];
plot_constellation(sig_mod(nsamples_per_symbol/2:nsamples_per_symbol:end),constellation_type,constellation_name,[-5:1:5]);
% Constellation



sig_mod_filt = ifft(fft(sig_mod).* fftshift(Hmf));
% Matched filtering.

eye = calc_eye(sig_mod_filt,nsamples_per_symbol);
% Calculate eye diagram.
figure('Name',['eye diagram (QAM16) after matched filter - in-phase: ' pulse_type])
plot(real(eye),'b','LineWidth',1.5);
xlabel('sample')
ylabel('amplitude (a.u.)')

figure('Name',['eye diagram (QAM16) after matched filter - quadrature: ' pulse_type])
plot(imag(eye),'b','LineWidth',1.5);
xlabel('sample')
ylabel('amplitude (a.u.)')

figure('Name',['eye diagram (QAM16) after matched filter - power: ' pulse_type])
plot(abs(eye).^2,'b','LineWidth',1.5);
xlabel('sample')
ylabel('amplitude (a.u.)')
% Plot eye diagram.



constellation_name = ['constellation (QAM16) after matched-filtering: ' pulse_type];
plot_constellation(sig_mod_filt(nsamples_per_symbol/2:nsamples_per_symbol:end),constellation_type,constellation_name,[-5:1:5]);
% Constellation


% We observe:
% 1. The pulses are properly timed, i.e. no delay is introduced by the
% matched filter.
% 2. The ISI that was introduced by the root raised-cosine pulse shaping is
% suppressed by the matched filter.


% -------------------------------------------------------------------------
% Check that the elec_filter_tf_matched is working fine
% -------------------------------------------------------------------------

% The matched filter calculation has been packaged into the
% elec_filter_tf_matched function.
% We check that it performs well.

tf = elec_tf_matched(sig_pulse,dt);
% Matched filter transfer function calculation.

sig_tf = elec_filter(sig_mod,tf);
% Matched filtering.

eye_tf = calc_eye(sig_tf,nsamples_per_symbol);
% Calculate eye diagram.
figure('Name',['eye diagram (QAM16) after functionalised matched filter - in-phase: ' pulse_type])
plot(real(eye_tf),'b','LineWidth',1.5);
xlabel('sample')
ylabel('amplitude (a.u.)')

figure('Name',['eye diagram (QAM16) after functionalised matched filter - quadrature: ' pulse_type])
plot(imag(eye_tf),'b','LineWidth',1.5);
xlabel('sample')
ylabel('amplitude (a.u.)')

figure('Name',['eye diagram (QAM16) after functionalised matched filter - power: ' pulse_type])
plot(abs(eye_tf).^2,'b','LineWidth',1.5);
xlabel('sample')
ylabel('amplitude (a.u.)')
% Plot eye diagram.



constellation_name = ['constellation (QAM16) after functionalised matched-filtering: ' pulse_type];
plot_constellation(sig_tf(nsamples_per_symbol/2:nsamples_per_symbol:end),constellation_type,constellation_name,[-5:1:5]);
% Constellation










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
