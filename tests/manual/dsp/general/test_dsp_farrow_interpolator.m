% -------------------------------------------------------------------------
% Test of polynomial interpolation using Farrow structure
% We consider the following interpolations:
% 1) (piecewise, Lagrange) linear (2 point basis)
% 2) (piecewise, Lagrange) parabolic (3 point basis)
% 3) (piecewise, Lagrange) cubic (4 point basis)
% 5) (piecewise) parabolic (4 point basis)
%
% We perform interpolation of a sinusoidal signal at a specified normalized
% frequency.
%
% The frequency responses (magnitude and phase delay) are also calculated.
% The presence of gain for the (piecewise) parabolic (3 point basis)
% interpolation for sufficiently large values of normalized frequencies can
% be observed.
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
time_stamp = datestr(datetime('now','TimeZone','Z'),'yyyymmddThhMMSSZ');

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
fig.interpreter = 'latex';


% -------------------------------------------------------------------------
% Test signal parameters
% -------------------------------------------------------------------------
f0 = 1;
% Signal frequency, in Hz.
phi = 0;
% Phase of the signal, in rad.


% -------------------------------------------------------------------------
% Time and frequency axes
% -------------------------------------------------------------------------
nsamples= 32;


% -------------------------------------------------------------------------
% Sampling
% -------------------------------------------------------------------------
fhat = 0.1
% Normalised frequency (cycles-per-sample).

nsamples_per_period = 1/fhat
% Number of samples per signal period (average).

fs = f0/fhat;
% Sampling frequency, in Hz.
Ts = 1/fs;
% Corresponding sampling period.


% -------------------------------------------------------------------------
% Generate signal samples
% -------------------------------------------------------------------------
sampling_instants = [0:1:nsamples - 1]*Ts;
x = cos(2*pi*f0*sampling_instants - phi);
% x(1:4) = 0;
% We can set the first samples of the signal to 0 to observe the delay, if
% needed.
time_array = linspace(0,sampling_instants(end),1000);
sig = cos(2*pi*f0*time_array - phi);

if do_debug    
    figure('Name','signal')
    plot(time_array,sig,'b-')
    hold on
    stem(sampling_instants,x,'sr-')
    xlabel('time')
    ylabel('amplitude')
    legend('original signal','samples')    
end


% -------------------------------------------------------------------------
% Fractional delay
% -------------------------------------------------------------------------
% mu = 0.5;
mu = rand(1,length(x));
% mu = 0.2*ones(1,length(x));
% mu(10) = 0.85;
% Fractional delay.



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Linear interpolator
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
fprintf('\n\n%s\n\n','Linear interpolator');

nfilters = 2;
ntaps = 2;
c0 = [0, 1];
c1 = [1, -1];
c = [c0;c1].';
% Coefficients of the Farrow structure for linear interpolation. 
 

y = dsp_farrow(x,c,zeros(ntaps,nfilters),mu);
% Interpolation with Farrow structure.

if do_debug
    figure('Name','linear - output samples')
    stem(x,'rs')
    hold on
    stem(y,'bo--')
    legend('input','output - Farrow')
    xlabel('sample index')
end

y = circshift(y,[0 -1]);
% Compensate delay, for representation.

figure('Name','linear - interpolated signal')
plot(time_array,sig,'b-')
hold on
stem(sampling_instants,x,'sr-')
stem(sampling_instants + mu*Ts,y,'ob--')
plot(sampling_instants,x,'r-')
xlabel('time')
ylabel('amplitude')
legend('original signal','samples','interpolated')


% -------------------------------------------------------------------------
% Frequency response of linear interpolator
% -------------------------------------------------------------------------
h00 = [0, 1];
h01 = [0.1, 1 - 0.1];
h02 = [0.2, 1 - 0.2];
h03 = [0.3, 1 - 0.3];
h04 = [0.4, 1 - 0.4];
h05 = [0.5, 1 - 0.5];
h06 = [0.6, 1 - 0.6];
h07 = [0.7, 1 - 0.7];
h08 = [0.8, 1 - 0.8];
h09 = [0.9, 1 - 0.9];
% Impulse responses for different values of mu within {0,0.9}

nfft = 512;
% Number of points for calculation of the filter response.
fnorm = [0:(nfft-1)/2]/nfft;
% Normalised frequencies.

H00 = fft(h00,nfft);
H01 = fft(h01,nfft);
H02 = fft(h02,nfft);
H03 = fft(h03,nfft);
H04 = fft(h04,nfft);
H05 = fft(h05,nfft);
H06 = fft(h06,nfft);
H07 = fft(h07,nfft);
H08 = fft(h08,nfft);
H09 = fft(h09,nfft);
% Frequency responses.

mag00 = 10*log10(abs(H00(1:nfft/2)).^2);
mag01 = 10*log10(abs(H01(1:nfft/2)).^2);
mag02 = 10*log10(abs(H02(1:nfft/2)).^2);
mag03 = 10*log10(abs(H03(1:nfft/2)).^2);
mag04 = 10*log10(abs(H04(1:nfft/2)).^2);
mag05 = 10*log10(abs(H05(1:nfft/2)).^2);
mag06 = 10*log10(abs(H06(1:nfft/2)).^2);
mag07 = 10*log10(abs(H07(1:nfft/2)).^2);
mag08 = 10*log10(abs(H08(1:nfft/2)).^2);
mag09 = 10*log10(abs(H09(1:nfft/2)).^2);
% Magnitude over positive frequencies,

phi00 = unwrap(angle(H00(1:nfft/2)));
phi01 = unwrap(angle(H01(1:nfft/2)));
phi02 = unwrap(angle(H02(1:nfft/2)));
phi03 = unwrap(angle(H03(1:nfft/2)));
phi04 = unwrap(angle(H04(1:nfft/2)));
phi05 = unwrap(angle(H05(1:nfft/2)));
phi06 = unwrap(angle(H06(1:nfft/2)));
phi07 = unwrap(angle(H07(1:nfft/2)));
phi08 = unwrap(angle(H08(1:nfft/2)));
phi09 = unwrap(angle(H09(1:nfft/2)));
% Phase over positive frequencies.

pd00 = -phi00./fnorm/2/pi;
pd01 = -phi01./fnorm/2/pi;
pd02 = -phi02./fnorm/2/pi;
pd03 = -phi03./fnorm/2/pi;
pd04 = -phi04./fnorm/2/pi;
pd05 = -phi05./fnorm/2/pi;
pd06 = -phi06./fnorm/2/pi;
pd07 = -phi07./fnorm/2/pi;
pd08 = -phi08./fnorm/2/pi;
pd09 = -phi09./fnorm/2/pi;
% Phase delay.


figure('Name','linear interpolator frequency response')
subplot(2,1,1)
plot(fnorm,mag00,'k-')
hold on
plot(fnorm,mag01,'b-')
plot(fnorm,mag02,'r-')
plot(fnorm,mag03,'g-')
plot(fnorm,mag04,'c-')
plot(fnorm,mag05,'m-')
plot(fnorm,mag06,'c--')
plot(fnorm,mag07,'g--')
plot(fnorm,mag08,'r--')
plot(fnorm,mag09,'b--')
ylim([-20 5])
xlabel('normalized frequency (cycles/sample)')
ylabel('magnitude (db)')
legend('\mu = 0','\mu = 0.1','\mu = 0.2','\mu = 0.3','\mu = 0.4','\mu = 0.5','\mu = 0.6','\mu = 0.7','\mu = 0.8','\mu = 0.9','Location','SouthWest')
subplot(2,1,2)
plot(fnorm,pd00,'k-')
hold on
plot(fnorm,pd01,'b-')
plot(fnorm,pd02,'r-')
plot(fnorm,pd03,'g-')
plot(fnorm,pd04,'c-')
plot(fnorm,pd05,'m-')
plot(fnorm,pd06,'c--')
plot(fnorm,pd07,'g--')
plot(fnorm,pd08,'r--')
plot(fnorm,pd09,'b--')
xlabel('normalized frequency (cycles/sample)')
ylabel('phase delay (sample)')
ylim([0 1.1])


%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Parabolic interpolator
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
fprintf('\n\n%s\n\n','Parabolic interpolator (3-basepoints): vorsicht, not linear phase!');

 nfilters = 3;
 ntaps = 3;
 c0 = [0,   1,  0];
 c1 = [0.5, 0,  -0.5];
 c2 = [0.5, -1, 0.5];

 c = [c0;c1;c2].';
 % Coefficients of the Farrow structure for parabolic interpolation. 
 
y = dsp_farrow(x,c,zeros(ntaps,nfilters),mu);
% Interpolation with Farrow structure.

if do_debug
    figure('Name','parabolic - output samples')
    stem(x,'rs')
    hold on
    stem(y,'bo--')
    legend('input','output - Farrow')
    xlabel('sample index')
end

y = circshift(y,[0 -1]);
% Compensate delay, for representation.

figure('Name','parabolic - interpolated signal')
plot(time_array,sig,'b-')
hold on
stem(sampling_instants,x,'sr-')
stem(sampling_instants + mu*Ts,y,'ob--')
xlabel('time')
ylabel('amplitude')
legend('original signal','samples','interpolated')


% -------------------------------------------------------------------------
% Response of parabolic interpolator
% -------------------------------------------------------------------------
mul = 0.0;
vmu = [1,mul,mul^2].';
h00 = c*vmu;
h00 = h00(:).';
% Filter coefficients for fractional delay mu = 0.

mul = 0.1;
vmu = [1,mul,mul^2].';
h01 = c*vmu;
h01 = h01(:).';
% Filter coefficients for fractional delay mu = 0.1.

mul = 0.2;
vmu = [1,mul,mul^2].';
h02 = c*vmu;
h02 = h02(:).';
% Filter coefficients for fractional delay mu = 0.2.

mul = 0.3;
vmu = [1,mul,mul^2].';
h03 = c*vmu;
h03 = h03(:).';
% Filter coefficients for fractional delay mu = 0.3.

mul = 0.4;
vmu = [1,mul,mul^2].';
h04 = c*vmu;
h04 = h04(:).';
% Filter coefficients for fractional delay mu = 0.4.

mul = 0.5;
vmu = [1,mul,mul^2].';
h05 = c*vmu;
h05 = h05(:).';
% Filter coefficients for fractional delay mu = 0.5.

mul = 0.6;
vmu = [1,mul,mul^2].';
h06 = c*vmu;
h06 = h06(:).';
% Filter coefficients for fractional delay mu = 0.6.

mul = 0.7;
vmu = [1,mul,mul^2].';
h07 = c*vmu;
h07 = h07(:).';
% Filter coefficients for fractional delay mu = 0.7.

mul = 0.8;
vmu = [1,mul,mul^2].';
h08 = c*vmu;
h08 = h08(:).';
% Filter coefficients for fractional delay mu = 0.8.

mul = 0.9;
vmu = [1,mul,mul^2].';
h09 = c*vmu;
h09 = h09(:).';
% Filter coefficients for fractional delay mu = 0.9.

nfft = 512;
% Number of points for calculation of the filter response.
fnorm = [0:(nfft-1)/2]/nfft;
% Normalised frequencies.

H00 = fft(h00,nfft);
H01 = fft(h01,nfft);
H02 = fft(h02,nfft);
H03 = fft(h03,nfft);
H04 = fft(h04,nfft);
H05 = fft(h05,nfft);
H06 = fft(h06,nfft);
H07 = fft(h07,nfft);
H08 = fft(h08,nfft);
H09 = fft(h09,nfft);
% Frequency responses.

mag00 = 10*log10(abs(H00(1:nfft/2)).^2);
mag01 = 10*log10(abs(H01(1:nfft/2)).^2);
mag02 = 10*log10(abs(H02(1:nfft/2)).^2);
mag03 = 10*log10(abs(H03(1:nfft/2)).^2);
mag04 = 10*log10(abs(H04(1:nfft/2)).^2);
mag05 = 10*log10(abs(H05(1:nfft/2)).^2);
mag06 = 10*log10(abs(H06(1:nfft/2)).^2);
mag07 = 10*log10(abs(H07(1:nfft/2)).^2);
mag08 = 10*log10(abs(H08(1:nfft/2)).^2);
mag09 = 10*log10(abs(H09(1:nfft/2)).^2);
% Magnitude over positive frequencies.

phi00 = unwrap(angle(H00(1:nfft/2)));
phi01 = unwrap(angle(H01(1:nfft/2)));
phi02 = unwrap(angle(H02(1:nfft/2)));
phi03 = unwrap(angle(H03(1:nfft/2)));
phi04 = unwrap(angle(H04(1:nfft/2)));
phi05 = unwrap(angle(H05(1:nfft/2)));
phi06 = unwrap(angle(H06(1:nfft/2)));
phi07 = unwrap(angle(H07(1:nfft/2)));
phi08 = unwrap(angle(H08(1:nfft/2)));
phi09 = unwrap(angle(H09(1:nfft/2)));
% Phase over positive frequencies.

pd00 = -phi00./fnorm/2/pi;
pd01 = -phi01./fnorm/2/pi;
pd02 = -phi02./fnorm/2/pi;
pd03 = -phi03./fnorm/2/pi;
pd04 = -phi04./fnorm/2/pi;
pd05 = -phi05./fnorm/2/pi;
pd06 = -phi06./fnorm/2/pi;
pd07 = -phi07./fnorm/2/pi;
pd08 = -phi08./fnorm/2/pi;
pd09 = -phi09./fnorm/2/pi;
% Phase delay.

figure('Name','parabolic interpolator frequency response')
subplot(2,1,1)
plot(fnorm,mag00,'k-')
hold on
plot(fnorm,mag01,'b-')
plot(fnorm,mag02,'r-')
plot(fnorm,mag03,'g-')
plot(fnorm,mag04,'c-')
plot(fnorm,mag05,'m-')
plot(fnorm,mag06,'c--')
plot(fnorm,mag07,'g--')
plot(fnorm,mag08,'r--')
plot(fnorm,mag09,'b--')
ylim([-20 5])
xlabel('normalized frequency (cycles/sample)')
ylabel('magnitude (db)')
legend('\mu = 0','\mu = 0.1','\mu = 0.2','\mu = 0.3','\mu = 0.4','\mu = 0.5','\mu = 0.6','\mu = 0.7','\mu = 0.8','\mu = 0.9','Location','SouthWest')
% subplot(3,1,2)
% plot(fnorm,phi00/pi,'k-')
% hold on
% plot(fnorm,phi01/pi,'b-')
% plot(fnorm,phi02/pi,'r-')
% plot(fnorm,phi03/pi,'g-')
% plot(fnorm,phi04/pi,'c-')
% plot(fnorm,phi05/pi,'m-')
% plot(fnorm,phi06/pi,'c--')
% plot(fnorm,phi07/pi,'g--')
% plot(fnorm,phi08/pi,'r--')
% plot(fnorm,phi09/pi,'b--')
% xlabel('normalized frequency (cycles/sample)')
% ylabel('phase (x\pi rad)')
subplot(2,1,2)
plot(fnorm,pd00,'k-')
hold on
plot(fnorm,pd01,'b-')
plot(fnorm,pd02,'r-')
plot(fnorm,pd03,'g-')
plot(fnorm,pd04,'c-')
plot(fnorm,pd05,'m-')
plot(fnorm,pd06,'c--')
plot(fnorm,pd07,'g--')
plot(fnorm,pd08,'r--')
plot(fnorm,pd09,'b--')
xlabel('normalized frequency (cycles/sample)')
ylabel('phase delay (sample)')
ylim([0 1.1])



%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Cubic Lagrange interpolator
% ------------------------------------------------------------------------- 
% ------------------------------------------------------------------------- 
fprintf('\n\n%s\n\n','Cubic interpolator');

nfilters = 4;
ntaps = 4;
c0 = [0,    0,     1,    0];
c1 = [-1/6, 1,     -1/2, -1/3];
c2 = [0,   1/2,    -1,   1/2];
c3 = [1/6, -1/2,   1/2,  -1/6];
 
c = [c0;c1;c2;c3].';
% Coefficients of the Farrow structure for cubic interpolation.
 
 
y = dsp_farrow(x,c,zeros(ntaps,nfilters),mu);
% Interpolation with Farrow structure.

if do_debug
    figure('Name','cubic - output samples')
    stem(x,'rs')
    hold on
    stem(y,'bo--')
    legend('input','output - Farrow')
    xlabel('sample index')
end

y = circshift(y,[0 -2]);
% Compensate delay, for representation.

figure('Name','cubic - interpolated signal')
plot(time_array,sig,'b-')
hold on
stem(sampling_instants,x,'sr-')
stem(sampling_instants + mu*Ts,y,'ob--')
xlabel('time')
ylabel('amplitude')
legend('original signal','samples','interpolated')

% -------------------------------------------------------------------------
% Response of cubic interpolator
% -------------------------------------------------------------------------
mul = 0.0;
vmu = [1,mul,mul^2,mul^3].';
h00 = c*vmu;
h00 = h00(:).';
% Filter coefficients for fractional delay mu = 0.

mul = 0.1;
vmu = [1,mul,mul^2,mul^3].';
h01 = c*vmu;
h01 = h01(:).';
% Filter coefficients for fractional delay mu = 0.1.

mul = 0.2;
vmu = [1,mul,mul^2,mul^3].';
h02 = c*vmu;
h02 = h02(:).';
% Filter coefficients for fractional delay mu = 0.2.

mul = 0.3;
vmu = [1,mul,mul^2,mul^3].';
h03 = c*vmu;
h03 = h03(:).';
% Filter coefficients for fractional delay mu = 0.3.

mul = 0.4;
vmu = [1,mul,mul^2,mul^3].';
h04 = c*vmu;
h04 = h04(:).';
% Filter coefficients for fractional delay mu = 0.4.

mul = 0.5;
vmu = [1,mul,mul^2,mul^3].';
h05 = c*vmu;
h05 = h05(:).';
% Filter coefficients for fractional delay mu = 0.5.

mul = 0.6;
vmu = [1,mul,mul^2,mul^3].';
h06 = c*vmu;
h06 = h06(:).';
% Filter coefficients for fractional delay mu = 0.6.

mul = 0.7;
vmu = [1,mul,mul^2,mul^3].';
h07 = c*vmu;
h07 = h07(:).';
% Filter coefficients for fractional delay mu = 0.7.

mul = 0.8;
vmu = [1,mul,mul^2,mul^3].';
h08 = c*vmu;
h08 = h08(:).';
% Filter coefficients for fractional delay mu = 0.8.

mul = 0.9;
vmu = [1,mul,mul^2,mul^3].';
h09 = c*vmu;
h09 = h09(:).';
% Filter coefficients for fractional delay mu = 0.9.

nfft = 512;
% Number of points for calculation of the filter response.
fnorm = [0:(nfft-1)/2]/nfft;
% Normalised frequencies.

H00 = fft(h00,nfft);
H01 = fft(h01,nfft);
H02 = fft(h02,nfft);
H03 = fft(h03,nfft);
H04 = fft(h04,nfft);
H05 = fft(h05,nfft);
H06 = fft(h06,nfft);
H07 = fft(h07,nfft);
H08 = fft(h08,nfft);
H09 = fft(h09,nfft);
% Frequency response.

mag00 = 10*log10(abs(H00(1:nfft/2)).^2);
mag01 = 10*log10(abs(H01(1:nfft/2)).^2);
mag02 = 10*log10(abs(H02(1:nfft/2)).^2);
mag03 = 10*log10(abs(H03(1:nfft/2)).^2);
mag04 = 10*log10(abs(H04(1:nfft/2)).^2);
mag05 = 10*log10(abs(H05(1:nfft/2)).^2);
mag06 = 10*log10(abs(H06(1:nfft/2)).^2);
mag07 = 10*log10(abs(H07(1:nfft/2)).^2);
mag08 = 10*log10(abs(H08(1:nfft/2)).^2);
mag09 = 10*log10(abs(H09(1:nfft/2)).^2);
% Magnitude for positive frequencies.

phi00 = unwrap(angle(H00(1:nfft/2)));
phi01 = unwrap(angle(H01(1:nfft/2)));
phi02 = unwrap(angle(H02(1:nfft/2)));
phi03 = unwrap(angle(H03(1:nfft/2)));
phi04 = unwrap(angle(H04(1:nfft/2)));
phi05 = unwrap(angle(H05(1:nfft/2)));
phi06 = unwrap(angle(H06(1:nfft/2)));
phi07 = unwrap(angle(H07(1:nfft/2)));
phi08 = unwrap(angle(H08(1:nfft/2)));
phi09 = unwrap(angle(H09(1:nfft/2)));
% Phase for positive frequencies.

pd00 = -phi00./fnorm/2/pi;
pd01 = -phi01./fnorm/2/pi;
pd02 = -phi02./fnorm/2/pi;
pd03 = -phi03./fnorm/2/pi;
pd04 = -phi04./fnorm/2/pi;
pd05 = -phi05./fnorm/2/pi;
pd06 = -phi06./fnorm/2/pi;
pd07 = -phi07./fnorm/2/pi;
pd08 = -phi08./fnorm/2/pi;
pd09 = -phi09./fnorm/2/pi;
% Phase delay.


figure('Name','cubic interpolator frequency response')
subplot(2,1,1)
plot(fnorm,mag00,'k-')
hold on
plot(fnorm,mag01,'b-')
plot(fnorm,mag02,'r-')
plot(fnorm,mag03,'g-')
plot(fnorm,mag04,'c-')
plot(fnorm,mag05,'m-')
plot(fnorm,mag06,'c--')
plot(fnorm,mag07,'g--')
plot(fnorm,mag08,'r--')
plot(fnorm,mag09,'b--')
ylim([-20 5])
xlabel('normalized frequency (cycles/sample)')
ylabel('magnitude (db)')
legend('\mu = 0','\mu = 0.1','\mu = 0.2','\mu = 0.3','\mu = 0.4','\mu = 0.5','\mu = 0.6','\mu = 0.7','\mu = 0.8','\mu = 0.9','Location','SouthWest')
subplot(2,1,2)
plot(fnorm,pd00,'k-')
hold on
plot(fnorm,pd01,'b-')
plot(fnorm,pd02,'r-')
plot(fnorm,pd03,'g-')
plot(fnorm,pd04,'c-')
plot(fnorm,pd05,'m-')
plot(fnorm,pd06,'c--')
plot(fnorm,pd07,'g--')
plot(fnorm,pd08,'r--')
plot(fnorm,pd09,'b--')
xlabel('normalized frequency (cycles/sample)')
ylabel('phase delay (sample)')
ylim([0.9 2.1])


%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Piecewise parabolic interpolator
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
fprintf('\n\n%s\n\n','Piecewise parabolic interpolator (4-basepoints)');

ntaps = 4;
nfilters = 3;
alpha = 1;
c = [0,  -alpha,     alpha;
     0,   alpha+1,  -alpha;
     1,   alpha-1   -alpha;
     0,  -alpha,     alpha];
% Coefficients of the Farrow structure for cubic interpolation.
 
 
y = dsp_farrow(x,c,zeros(ntaps,nfilters),mu);
% Interpolation with Farrow structure.

if do_debug
    figure('Name','piecewise-parabolic - output samples')
    stem(x,'rs')
    hold on
    stem(y,'bo--')
    legend('input','output - Farrow')
    xlabel('sample index')
end

y = circshift(y,[0 -2]);
% Compensate delay, for representation.

figure('Name','piecewise-parabolic - interpolated signal')
plot(time_array,sig,'b-')
hold on
stem(sampling_instants,x,'sr-')
stem(sampling_instants + mu*Ts,y,'ob--')
xlabel('time')
ylabel('amplitude')
legend('original signal','samples','interpolated')


% -------------------------------------------------------------------------
% Response of piecewise parabolic interpolator
% -------------------------------------------------------------------------
mul = 0.0;
vmu = [1,mul,mul^2].';
h00 = c*vmu;
h00 = h00(:).';
% Filter coefficients for fractional delay mu = 0.

mul = 0.1;
vmu = [1,mul,mul^2].';
h01 = c*vmu;
h01 = h01(:).';
% Filter coefficients for fractional delay mu = 0.1.

mul = 0.2;
vmu = [1,mul,mul^2].';
h02 = c*vmu;
h02 = h02(:).';
% Filter coefficients for fractional delay mu = 0.2.

mul = 0.3;
vmu = [1,mul,mul^2].';
h03 = c*vmu;
h03 = h03(:).';
% Filter coefficients for fractional delay mu = 0.3.

mul = 0.4;
vmu = [1,mul,mul^2].';
h04 = c*vmu;
h04 = h04(:).';
% Filter coefficients for fractional delay mu = 0.4.

mul = 0.5;
vmu = [1,mul,mul^2].';
h05 = c*vmu;
h05 = h05(:).';
% Filter coefficients for fractional delay mu = 0.5.

mul = 0.6;
vmu = [1,mul,mul^2].';
h06 = c*vmu;
h06 = h06(:).';
% Filter coefficients for fractional delay mu = 0.6.

mul = 0.7;
vmu = [1,mul,mul^2].';
h07 = c*vmu;
h07 = h07(:).';
% Filter coefficients for fractional delay mu = 0.7.

mul = 0.8;
vmu = [1,mul,mul^2].';
h08 = c*vmu;
h08 = h08(:).';
% Filter coefficients for fractional delay mu = 0.8.

mul = 0.9;
vmu = [1,mul,mul^2].';
h09 = c*vmu;
h09 = h09(:).';
% Filter coefficients for fractional delay mu = 0.9.

nfft = 512;
% Number of points for calculation of the filter response.
fnorm = [0:(nfft-1)/2]/nfft;
% Normalised frequencies.

H00 = fft(h00,nfft);
H01 = fft(h01,nfft);
H02 = fft(h02,nfft);
H03 = fft(h03,nfft);
H04 = fft(h04,nfft);
H05 = fft(h05,nfft);
H06 = fft(h06,nfft);
H07 = fft(h07,nfft);
H08 = fft(h08,nfft);
H09 = fft(h09,nfft);
% Frequency response.

mag00 = 10*log10(abs(H00(1:nfft/2)).^2);
mag01 = 10*log10(abs(H01(1:nfft/2)).^2);
mag02 = 10*log10(abs(H02(1:nfft/2)).^2);
mag03 = 10*log10(abs(H03(1:nfft/2)).^2);
mag04 = 10*log10(abs(H04(1:nfft/2)).^2);
mag05 = 10*log10(abs(H05(1:nfft/2)).^2);
mag06 = 10*log10(abs(H06(1:nfft/2)).^2);
mag07 = 10*log10(abs(H07(1:nfft/2)).^2);
mag08 = 10*log10(abs(H08(1:nfft/2)).^2);
mag09 = 10*log10(abs(H09(1:nfft/2)).^2);
% Magnitude for positive frequencies.

phi00 = unwrap(angle(H00(1:nfft/2)));
phi01 = unwrap(angle(H01(1:nfft/2)));
phi02 = unwrap(angle(H02(1:nfft/2)));
phi03 = unwrap(angle(H03(1:nfft/2)));
phi04 = unwrap(angle(H04(1:nfft/2)));
phi05 = unwrap(angle(H05(1:nfft/2)));
phi06 = unwrap(angle(H06(1:nfft/2)));
phi07 = unwrap(angle(H07(1:nfft/2)));
phi08 = unwrap(angle(H08(1:nfft/2)));
phi09 = unwrap(angle(H09(1:nfft/2)));
% Phase for positive frequencies.

pd00 = -phi00./fnorm/2/pi;
pd01 = -phi01./fnorm/2/pi;
pd02 = -phi02./fnorm/2/pi;
pd03 = -phi03./fnorm/2/pi;
pd04 = -phi04./fnorm/2/pi;
pd05 = -phi05./fnorm/2/pi;
pd06 = -phi06./fnorm/2/pi;
pd07 = -phi07./fnorm/2/pi;
pd08 = -phi08./fnorm/2/pi;
pd09 = -phi09./fnorm/2/pi;
% Phase delay.


figure('Name','piecewise parabolic interpolator frequency response')
subplot(2,1,1)
plot(fnorm,mag00,'k-')
hold on
plot(fnorm,mag01,'b-')
plot(fnorm,mag02,'r-')
plot(fnorm,mag03,'g-')
plot(fnorm,mag04,'c-')
plot(fnorm,mag05,'m-')
plot(fnorm,mag06,'c--')
plot(fnorm,mag07,'g--')
plot(fnorm,mag08,'r--')
plot(fnorm,mag09,'b--')
ylim([-20 5])
xlabel('normalized frequency (cycles/sample)')
ylabel('magnitude (db)')
legend('\mu = 0','\mu = 0.1','\mu = 0.2','\mu = 0.3','\mu = 0.4','\mu = 0.5','\mu = 0.6','\mu = 0.7','\mu = 0.8','\mu = 0.9','Location','SouthWest')
subplot(2,1,2)
plot(fnorm,pd00,'k-')
hold on
plot(fnorm,pd01,'b-')
plot(fnorm,pd02,'r-')
plot(fnorm,pd03,'g-')
plot(fnorm,pd04,'c-')
plot(fnorm,pd05,'m-')
plot(fnorm,pd06,'c--')
plot(fnorm,pd07,'g--')
plot(fnorm,pd08,'r--')
plot(fnorm,pd09,'b--')
xlabel('normalized frequency (cycles/sample)')
ylabel('phase delay (sample)')
ylim([0.9 2.1])


% ------
% Samples for fixed values of mu
% ------
if do_debug
    
    mu_fixed = 0;
    y = dsp_farrow(x,c,zeros(ntaps,nfilters),mu_fixed);
    
    figure('Name',['piecewise-parabolic - output samples for mu = ' num2str(mu_fixed)])
    stem(x,'rs')
    hold on
    stem(y,'bo--')
    legend('input','output')
    title(['\mu = ' num2str(mu_fixed)])
    % mu_fixed = 0
    % The interpolator does nothing.
    % We just observe a 2-sample delay.
    
    mu_fixed = 0.5;
    y = dsp_farrow(x,c,zeros(ntaps,nfilters),mu_fixed);
    
    figure('Name',['piecewise-parabolic - output samples for mu = ' num2str(mu_fixed)])
    stem(x,'rs')
    hold on
    stem(y,'bo--')
    legend('input','output')
    title(['\mu = ' num2str(mu_fixed)])
    
end



%%
% -------------------------------------------------------------------------
% Continuous time impulse responses of linear and cubic interpolators
% -------------------------------------------------------------------------
if do_debug
    
    nsa = 1024;
    % Number of points.
    
    mmax = 15;
    % Range.
    
    tn = linspace(-mmax ,mmax ,nsa);
    % Normalized time axis.
    
    fn = (0:nsa-1)/2/mmax;
    %  Normalized frequency axis.
    
    ha_lin = zeros(1,length(tn));    
    ha_lin( tn >= -1) = 1 + tn(tn >= -1);
    ha_lin( tn >= 0) = 1 - tn(tn >= 0);
    ha_lin( tn >= 1) = 0;
    % Linear response.
    
    
    ha_cub = zeros(1,length(tn));    
    ha_cub( tn >= -2) = (1/6)*(2 + tn(tn >= -2)).^3 - (1/6)*(2 + tn(tn >= -2));
    ha_cub( tn >= -1) = -(1/2)*(1 + tn(tn >= -1)).^3 + (1/2)*(1 + tn(tn >= -1)).^2 + (1 + tn(tn >= -1));
    ha_cub( tn >= 0) = (1/2)*(0 + tn(tn >= 0)).^3 - (0 + tn(tn >= 0)).^2 -(1/2)*(0 + tn(tn >= 0)) + 1;
    ha_cub( tn >= 1) = -(1/6)*(-1 + tn(tn >= 1)).^3 + (1/2)*(-1 + tn(tn >= 1)).^2 -(1/3)*(-1 + tn(tn >= 1));
    ha_cub( tn >= 2) = 0;
    % Cubic response.
    
    figure('Name','impulse response of linear and cubic Lagrange interpolators')
    plot(tn,ha_lin,'b-')
    hold on
    plot(tn,ha_cub,'r-')
    xlabel('t/T_{in}')
    ylabel('h_a(t)')
    legend('linear','cubic')
    xlim([-3 3])
    
    
    Ha_cub = fft(ha_cub)/length(ha_cub);
    Ha_lin = fft(ha_lin)/length(ha_lin);
    % Frequency responses.
    
    figure('Name','frequency response of linear and cubic Lagrange interpolators')
    plot(fn,10*log10(abs(Ha_lin).^2/max(abs(Ha_lin).^2)),'b-')
    hold on
    plot(fn,10*log10(abs(Ha_cub).^2/max(abs(Ha_cub).^2)),'r-')
    xlabel('f T_{in}')
    ylabel('H_a(f) (dB)')
    legend('linear','cubic')
    xlim([0,2.5])
    ylim([-60 5]);
    
end

% -------------------------------------------------------------------------
% Figure alignment
% -------------------------------------------------------------------------
% Requires alignfigs function
% Available at https://github.com/nickhale/alignfigs
alignfigs(3)

% -------------------------------------------------------------------------
% Save data
% -------------------------------------------------------------------------
% file_name_data = [file_name_core_data '_' time_stamp]
% save(file_name_data,'');

