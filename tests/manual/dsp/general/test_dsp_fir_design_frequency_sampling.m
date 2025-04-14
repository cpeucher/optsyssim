% -------------------------------------------------------------------------
% FIR filter design using the frequency sampling method
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Clean-up
% -------------------------------------------------------------------------
clear all
close all
format long

% -------------------------------------------------------------------------
% Strings for filenames
% -------------------------------------------------------------------------
file_name_core_figure = strrep(mfilename,'run','fig');
file_name_core_data = strrep(mfilename,'run','data');
time_stamp = string(datetime('now','TimeZone','Z','Format','yyyyMMdd''T''HHmmss''Z'));

% ------------------------------------------------------------------------- 
% Reinitialise the random number generator for reproducibility.
% ------------------------------------------------------------------------- 
stream = RandStream.getGlobalStream;
reset(stream);

% -------------------------------------------------------------------------
% Switches
% -------------------------------------------------------------------------
do_debug = 0;
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
% Parameters
% -------------------------------------------------------------------------
fs = 120e9;
% Sampling frequency, in Hz. This is chosen arbitrarily.

freq = linspace(0,fs/2,1000);
% Frequency axis
fnorm = freq/fs;
% Corresponding normalized frequency axis



% -------------------------------------------------------------------------
% Test 1: Rectangular filter with cutoff at pi/4 (normalized angular
% frequency) - 9 tap FIR filter
% -------------------------------------------------------------------------
fc = 0.25*pi/(2*pi)*fs;
% Cut-off frequency. Corresponds to pi/4 in normalized angular frequency.

params_elpf.type = 'rectangular';
params_elpf.f3dB = fc;
tf = elec_tf_elpf(params_elpf,freq);
% Calculate analog transfer function.

if do_debug
    figure('Name','analog transfer function')
    plot(freq/1.0e9,abs(tf).^2)
    xlabel('frequency (GHz)')
    ylabel('|H|^2')
end


N = 9;
% Number of taps. Choose N odd.

fsamp_norm = [0:1:(N-1)/2]/N;
% Normalized sampling frequencies over [0,1/2] interval.
fsamp = fsamp_norm*fs; 
% Corresponding frequencies.

Hsamp = elec_tf_elpf(params_elpf,fsamp);
% Acquire samples over [0,1/2] normalized frequency interval.

Hsamp = [Hsamp,fliplr(conj(Hsamp(2:end)))];
% Reconstitute samples of full frequency response over [0,1[ normalized frequency
% interval.

h = ifft(Hsamp);
h = fftshift(h);
% Filter tap coefficients.

figure('Name','tap coefficients - rectangular - 9 taps')
stem([0:N-1],h)
xlabel('n')

wh = [2*pi*fnorm];
% Normalized angular frequency over [0,pi] interval.
kk = [0:1:N-1].';
H = h*exp(-1i*kk*wh);
% Calculate the frequency response over [0,1/2] normalized frequency
% interval.

figure('Name','analog and digital filter responses - rectangular - 9 taps')
plot(2*fnorm,abs(tf).^2,'b-')
hold on
plot(2*fsamp_norm,abs(Hsamp(1:1+(N-1)/2)).^2,'ro')
plot(2*fnorm,abs(H).^2,'r-')
xlabel('normalized angular frequency (x\pi radians-per-sample)')
ylabel('|H|^2')
legend('desired analog filter response','samples','digital filter response')


% We created a frequency sampling filter design function on 2021-10-05.
% We test it below
h2 = dsp_fir_design_frequency_sampling(freq,tf,fs,N,1);
% OK.


%%
% -------------------------------------------------------------------------
% Test 2: Rectangular filter with cutoff at pi/4 (normalized angular
% frequency) - 25 tap FIR filter
% -------------------------------------------------------------------------
fc = 0.25*pi/(2*pi)*fs;
% Cut-off frequency. Corresponds to pi/4 in normalized angular frequency.

params_elpf.type = 'rectangular';
params_elpf.f3dB = fc;
tf = elec_tf_elpf(params_elpf,freq);
% Calculate analog transfer function.

if do_debug
    figure('Name','analog transfer function')
    plot(freq/1.0e9,abs(tf).^2)
    xlabel('frequency (GHz)')
    ylabel('|H|^2')
end


N = 25;
% Number of taps. Choose N odd.

fsamp_norm = [0:1:(N-1)/2]/N;
% Normalized sampling frequencies over [0,1/2] interval.
fsamp = fsamp_norm*fs; 
% Corresponding frequencies.

Hsamp = elec_tf_elpf(params_elpf,fsamp);
% Acquire samples over [0,1/2] normalized frequency interval.

Hsamp = [Hsamp,fliplr(conj(Hsamp(2:end)))];
% Reconstitute samples of full frequency response over [0,1[ normalized frequency
% interval.

h = ifft(Hsamp);
h = fftshift(h);
% Filter tap coefficients.

figure('Name','tap coefficients - - rectangular - 25 taps')
stem([0:N-1],h)
xlabel('n')

wh = [2*pi*fnorm];
% Normalized angular frequency over [0,pi] interval.
kk = [0:1:N-1].';
H = h*exp(-1i*kk*wh);
% Calculate the frequency response over [0,1/2] normalized frequency
% interval.

figure('Name','analog and digital filter responses - rectangular - 25 taps')
plot(2*fnorm,abs(tf),'b-')
hold on
plot(2*fsamp_norm,abs(Hsamp(1:1+(N-1)/2)),'ro')
plot(2*fnorm,abs(H),'r-')
xlabel('normalized angular frequency (x\pi radians-per-sample)')
ylabel('|H|')
legend('desired analog filter response','samples','digital filter response')


%%
% -------------------------------------------------------------------------
% Test 3: Butterworth filter
% -------------------------------------------------------------------------
params_elpf.type = 'butterworth';%'bessel';
params_elpf.order = 4;
params_elpf.f3dB = fs/4;
tf = elec_tf_elpf(params_elpf,freq);
% Calculate analog transfer function.

tf_phase = unwrap(angle(tf));
tf_delay_phase =  -tf_phase./freq/2/pi;
[freq2,tf_delay_group] = num_diff(1,freq,tf_phase);
tf_delay_group = -tf_delay_group/2/pi;
% Calculate phase and group delay of the filter.



if do_debug
    figure('Name','analog transfer function')
    plot(freq/1.0e9,abs(tf).^2)
    xlabel('frequency (GHz)')
    ylabel('|H|^2')
end


N = 25;
% Number of taps. Choose N odd.

fsamp_norm = [0:1:(N-1)/2]/N;
% Normalized sampling frequencies over [0,1/2] interval.
fsamp = fsamp_norm*fs; 
% Corresponding frequencies.

Hsamp = elec_tf_elpf(params_elpf,fsamp);
% Acquire samples over [0,1/2] normalized frequency interval.

Hsamp = [Hsamp,fliplr(conj(Hsamp(2:end)))];
% Reconstitute samples of full frequency response over [0,1[ normalized frequency
% interval.

h = ifft(Hsamp);
h = fftshift(h);
% Filter tap coefficients.

figure('Name','tap coefficients - butterworth - 25 taps')
stem([0:N-1],h)
xlabel('n')

wh = [2*pi*fnorm];
% Normalized angular frequency over [0,pi] interval.
kk = [0:1:N-1].';
H = h*exp(-1i*kk*wh);
% Calculate the frequency response over [0,1/2] normalized frequency
% interval.



H_phase = unwrap(angle(H));
H_delay_phase =  -H_phase./fnorm/2/pi;
[fnorm2,H_delay_group] = num_diff(1,fnorm,H_phase);
H_delay_group = -H_delay_group/2/pi;
% Calculate phase and group delay of the filter.

figure('Name','analog and digital filter responses - butterworth - 25 taps')
subplot(2,1,1)
plot(2*fnorm,abs(tf),'b-')
hold on
plot(2*fsamp_norm,abs(Hsamp(1:1+(N-1)/2)),'ro')
plot(2*fnorm,abs(H),'r-')
xlabel('normalized angular frequency (x\pi radians-per-sample)')
ylabel('|H|')
legend('desired analog filter response','samples','digital filter response')
subplot(2,1,2)
plot(2*fnorm2,tf_delay_group,'b-')
hold on
plot(2*fnorm2,H_delay_group/fs-(N-1)/2/fs,'r-')
xlabel('normalized angular frequency (x\pi radians-per-sample)')
ylabel('relative group delay (s)')


%%
% -------------------------------------------------------------------------
% Test 4: Bessel filter
% -------------------------------------------------------------------------
params_elpf.type = 'bessel';
params_elpf.order = 4;
params_elpf.f3dB = fs/4;
tf = elec_tf_elpf(params_elpf,freq);
% Calculate analog transfer function.

tf_phase = unwrap(angle(tf));
tf_delay_phase =  -tf_phase./freq/2/pi;
[freq2,tf_delay_group] = num_diff(1,freq,tf_phase);
tf_delay_group = -tf_delay_group/2/pi;
% Calculate phase and group delay of the filter.



if do_debug
    figure('Name','analog transfer function')
    plot(freq/1.0e9,abs(tf).^2)
    xlabel('frequency (GHz)')
    ylabel('|H|^2')
end


N = 25;
% Number of taps. Choose N odd.

fsamp_norm = [0:1:(N-1)/2]/N;
% Normalized sampling frequencies over [0,1/2] interval.
fsamp = fsamp_norm*fs; 
% Corresponding frequencies.

Hsamp = elec_tf_elpf(params_elpf,fsamp);
% Acquire samples over [0,1/2] normalized frequency interval.

Hsamp = [Hsamp,fliplr(conj(Hsamp(2:end)))];
% Reconstitute samples of full frequency response over [0,1[ normalized frequency
% interval.

h = ifft(Hsamp);
h = fftshift(h);
% Filter tap coefficients.

figure('Name','tap coefficients - bessel - 25 taps')
stem([0:N-1],h)
xlabel('n')

wh = [2*pi*fnorm];
% Normalized angular frequency over [0,pi] interval.
kk = [0:1:N-1].';
H = h*exp(-1i*kk*wh);
% Calculate the frequency response over [0,1/2] normalized frequency
% interval.



H_phase = unwrap(angle(H));
H_delay_phase =  -H_phase./fnorm/2/pi;
[fnorm2,H_delay_group] = num_diff(1,fnorm,H_phase);
H_delay_group = -H_delay_group/2/pi;
% Calculate phase and group delay of the filter.

figure('Name','analog and digital filter responses - bessel - 25 taps')
subplot(2,1,1)
plot(2*fnorm,abs(tf),'b-')
hold on
plot(2*fsamp_norm,abs(Hsamp(1:1+(N-1)/2)),'ro')
plot(2*fnorm,abs(H),'r-')
xlabel('normalized angular frequency (x\pi radians-per-sample)')
ylabel('|H|')
legend('desired analog filter response','samples','digital filter response')
subplot(2,1,2)
plot(2*fnorm2,tf_delay_group,'b-')
hold on
plot(2*fnorm2,H_delay_group/fs-(N-1)/2/fs,'r-')
xlabel('normalized angular frequency (x\pi radians-per-sample)')
ylabel('relative group delay (s)')



%%
% -------------------------------------------------------------------------
% Test 5: digital filter for pre-equalization of the response of zero-order
% hold circuit in a DAC.
% -------------------------------------------------------------------------

tf = 1./func_sinc(pi*freq/fs);
% Required compensation transfer function


figure('Name','analog transfer function')
plot(2*fnorm,10*log10(abs(tf).^2),'b-')
xlabel('normalized angular frequency (x\pi radians-per-sample)')
ylabel('|H|^2 (dB)')

N = 25;
% Number of taps. Choose N odd.

fsamp_norm = [0:1:(N-1)/2]/N;
% Normalized sampling frequencies over [0,1/2] interval.
fsamp = fsamp_norm*fs; 
% Corresponding frequencies.

Hsamp = 1./func_sinc(pi*fsamp/fs);
% Acquire samples over [0,1/2] normalized frequency interval.

Hsamp = [Hsamp,fliplr(conj(Hsamp(2:end)))];
% Reconstitute samples of full frequency response over [0,1[ normalized frequency
% interval.

h = ifft(Hsamp);
h = fftshift(h);
% Filter tap coefficients.

figure('Name','tap coefficients - ZOH compensation - 25 taps')
stem([0:N-1],h)
xlabel('n')

wh = [2*pi*fnorm];
% Normalized angular frequency over [0,pi] interval.
kk = [0:1:N-1].';
H = h*exp(-1i*kk*wh);
% Calculate the frequency response over [0,1/2] normalized frequency
% interval.


figure('Name','analog and digital filter responses - ZOH compensation - 25 taps')
plot(2*fnorm,10*log10(abs(tf).^2),'b-')
hold on
plot(2*fsamp_norm,10*log10(abs(Hsamp(1:1+(N-1)/2)).^2),'ro')
plot(2*fnorm,10*log10(abs(H).^2),'r-')
xlabel('normalized angular frequency (x\pi radians-per-sample)')
ylabel('|H|^2 (dB)')
legend('desired analog filter response','samples','digital filter response')




%%
% -------------------------------------------------------------------------
% Test 6: Example from 
% J.G. Proakis and D.G. Manolakis, Digital Signal Processing, 
% Fourth edition, Pearson new international edition 
%(Pearson, Harlow, Essex, 2014).
% Example 2.1, page 688, chapter 10.
% -------------------------------------------------------------------------

N = 15;
% Number of taps. Choose N odd.

fsamp_norm = [0:1:(N-1)/2]/N;
% Normalized sampling frequencies over [0,1/2] interval.
fsamp = fsamp_norm*fs; 
% Corresponding frequencies.

Hsamp = [1 1 1 1 0.4 0 0 0];
% Samples over [0,1/2] normalized frequency interval.

Hsamp = [Hsamp,fliplr(conj(Hsamp(2:end)))];
% Reconstitute samples of full frequency response over [0,1[ normalized frequency
% interval.

h = ifft(Hsamp);
h = fftshift(h)
% Filter tap coefficients.

figure('Name','tap coefficients - Proakis 2.1')
stem([0:N-1],h)
xlabel('n')

wh = [2*pi*fnorm];
% Normalized angular frequency over [0,pi] interval.
kk = [0:1:N-1].';
H = h*exp(-1i*kk*wh);
% Calculate the frequency response over [0,1/2] normalized frequency
% interval.



figure('Name','analog and digital filter responses - Proakis 2.1')
plot(2*fsamp_norm,abs(Hsamp(1:1+(N-1)/2)).^2,'ro')
hold on
plot(2*fnorm,abs(H).^2,'r-')
xlabel('normalized angular frequency (x\pi radians-per-sample)')
ylabel('|H|^2')
legend('samples','digital filter response')

figure('Name','analog and digital filter responses - Proakis 2.1 - log')
plot(2*fsamp_norm,10*log10(abs(Hsamp(1:1+(N-1)/2)).^2),'ro')
hold on
plot(2*fnorm,10*log10(abs(H).^2),'r-')
xlabel('normalized angular frequency (x\pi radians-per-sample)')
ylabel('|H|^2 (dB)')
legend('samples','digital filter response')
ylim([-80 10]);


%%
% -------------------------------------------------------------------------
% Test 7: We add a window function to the previous example.
% -------------------------------------------------------------------------

params_window.type = 'kaiser';%'hamming';'hann';'bartlett';'rectangular';
params_window.length = 15;
params_window.symmetry = 'periodic';%'symmetric';
params_window.beta = 2.4*pi;% For Kaiser window.
w = dsp_window(params_window);

h = w.*h;
% Apply window.

[~,H] = dsp_fir_frequency_response(h);
% Calculate frequency response.

figure('Name','analog and digital filter responses - Proakis 2.1 - windowed - log')
plot(2*fsamp_norm,10*log10(abs(Hsamp(1:1+(N-1)/2)).^2),'ro')
hold on
plot(2*fnorm,10*log10(abs(H).^2),'r-')
xlabel('normalized angular frequency (x\pi radians-per-sample)')
ylabel('|H|^2 (dB)')
legend('samples','digital filter response')
ylim([-80 10]);


%%
% -------------------------------------------------------------------------
% Figure alignment
% -------------------------------------------------------------------------
% Requires alignfigs function
% Available at https://github.com/nickhale/alignfigs
alignfigs(5)


% -------------------------------------------------------------------------
% Save data
% -------------------------------------------------------------------------
% file_name_data = [file_name_core_data '_' time_stamp]
% save(file_name_data,'');