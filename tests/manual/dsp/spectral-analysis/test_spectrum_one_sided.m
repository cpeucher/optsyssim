% -------------------------------------------------------------------------
% Test of the spectrum_one_sided.m function.
%
% 2021-05-19
%
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Clean-up
% -------------------------------------------------------------------------
clear all
close all


% -------------------------------------------------------------------------
% String for filenames
% -------------------------------------------------------------------------
file_name_figure_core = strrep(mfilename,'run','fig');
file_name_data_core = strrep(mfilename,'run','data');
time_stamp = datestr(datetime('now','TimeZone','Z'),'yyyymmddThhMMSSZ');
% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% Global parameters
% -------------------------------------------------------------------------
nsamples = 2^10;
% Number of samples of the signal.

fs = 300;
% Sampling frequency

dt = 1/fs;
% Interval between consecutive samples in the time domain.
df = fs/nsamples;
% Interval between consecutive samples in the frequency domain.

time_array = (0:nsamples-1)*dt;
% Construct time vector.
frequency_array = (-nsamples/2:nsamples/2-1)*df;
% Construct frequency vector.


% -------------------------------------------------------------------------
% Test signal
% -------------------------------------------------------------------------
xdc = 4;
% DC level of the signal.
xac1 = 2;
% Amplitude of the cosine component.
xac2 = 0.7;
% Amplitude of the sine component.
f1 = floor(3/df)*df
% Frequency of the cosine component. 
% We want to make sure it is on the grid.
f2 = floor(4/df)*df
% Frequency of the sine component.
% We want to make sure it is on the grid.
dphi = 3*pi/8
% Relative phase between the two terms.


x = xdc + xac1*cos(2*pi*f1*time_array) + xac2*cos(2*pi*f2*time_array - dphi);
% Test signal.

figure('Name','test signal')
plot(time_array,x)
xlabel('time (s)')
ylabel('voltage (a.u.)')


% -------------------------------------------------------------------------
% Power level of the different frequency components
% -------------------------------------------------------------------------
dc_power_expected_db = 10*log10(xdc^2);
ac1_power_expected_db = 10*log10(xac1^2/2);
ac2_power_expected_db = 10*log10(xac2^2/2);


% -------------------------------------------------------------------------
% Double-sided spectrum
% -------------------------------------------------------------------------

X2 = fftshift(fft(x)/nsamples);

figure('Name','double-sided power spectrum')
plot(frequency_array,10*log10(abs(X2).^2))
xlabel('frequency (Hz)')
ylabel('power (dB)')
xlim([-5 5])
hline(dc_power_expected_db,'k--',['DC power = ' num2str(dc_power_expected_db) ' dB'])
hline(ac1_power_expected_db - 10*log10(2),'r--',['AC1 power - 3dB = ' num2str(ac1_power_expected_db - 10*log10(2)) ' dB'])
hline(ac2_power_expected_db - 10*log10(2),'b--',['AC2 power - 3dB = ' num2str(ac2_power_expected_db - 10*log10(2)) ' dB'])


% -------------------------------------------------------------------------
% Single sided spectrum with nfft = nsamples
% -------------------------------------------------------------------------
[f,X1] = spectrum_one_sided(x,fs);

figure('Name','double-sided power spectrum with nfft = nsamples')
plot(f,10*log10(abs(X1).^2))
xlabel('frequency (Hz)')
ylabel('power (dB)')
xlim([0 5])
hline(dc_power_expected_db,'k--',['DC power = ' num2str(dc_power_expected_db) ' dB'])
hline(ac1_power_expected_db,'r--',['AC1 power= ' num2str(ac1_power_expected_db) ' dB'])
hline(ac2_power_expected_db,'b--',['AC2 power = ' num2str(ac2_power_expected_db) ' dB'])



X1f1 = X1(find(f == f1));
X1f2 = X1(find(f == f2));
dphi_retrieved = angle(X1f1) - angle(X1f2)
% Check that we can retrieve the phase.




% -------------------------------------------------------------------------
% Save data
% -------------------------------------------------------------------------
% file_name_data = [file_name_data_core '_' time_stamp]
% save(file_name_data,'');




