% -------------------------------------------------------------------------
% We test the implementation of the standard carrier frequency offset 
% estimation algorithm described in
% A. Leven, N. Kaneda, U.-V. Koc, and Y.-K. Chen, "Frequency estimation in 
% intradyne reception," IEEE Photonics Technology Letters 19, 
% 366?368 (2007) [DOI: 10.1109/LPT.2007.891893].
%
% We generate a (discrete-time) QPSK signal, add some frequency offset, 
% possibly some AWGN, then apply the cfo_estimation_leven.m function to
% estimate the CFO.
%
% The added CFO can be constant or sinusoidal.
% 
% Estimation is performed over 2 samples separated by a specified delay 
% (point) or averaged over a block of 2 samples separated by a specified 
% delay (block).
%
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Clean-up
% -------------------------------------------------------------------------
clear all
close all

% -------------------------------------------------------------------------
% Reset random numbers generator
% -------------------------------------------------------------------------
stream = RandStream.getGlobalStream;
reset(stream);

% -------------------------------------------------------------------------
% Global parameters
% -------------------------------------------------------------------------
nsymbols = 2^17;
% Number of symbols. Should be power of 2.

% -------------------------------------------------------------------------
% Define constellation
% -------------------------------------------------------------------------
m = 4;
% Constellation order. M points with log2(M) even (square constellation)

[constellation,norm_es,norm_emax] = define_constellation('qpsk_gray',m);
% Define the constellation for Gray mapping

% -------------------------------------------------------------------------
% Binary data stream
% -------------------------------------------------------------------------
data = generate_binary(nsymbols,m);
% Generate binary data

[words_dec,words_bin] = conv_bin2dec(data,log2(m));
% Convert binary data into 2-bit words

% -------------------------------------------------------------------------
% Mapping
% -------------------------------------------------------------------------
symbs = mapping(words_dec,constellation);
% Mapping

symbs_ref = symbs;
% Save reference signal

% -------------------------------------------------------------------------
% Check constellation at transmitter output
% -------------------------------------------------------------------------
plot_constellation(symbs,'plain','Constellation at Tx output');

% -------------------------------------------------------------------------
% Parameters
% ------------------------------------------------------------------------- 
symbol_rate = 25e9;     
% Symbol rate, in baud. Signal is sampled at the symbol rate.
esn0_db = +Inf;
% esn0_db = 30;
% Signal-to-noise ratio, in dB. Will govern the amount of AWGN added.

% -------------------------------------------------------------------------
% Added CFO parameters
% -------------------------------------------------------------------------
cfo_mean = 300e6;       % Mean value of the CFO, in Hz
cfo_pp = 1*10e6;        % Peak-to-peak value of the sinusoidal CFO
cfo_var_freq = symbol_rate/20000; % Frequency of the variations of the CFO, when sinusoidal

% -------------------------------------------------------------------------
% CFO estimation parameters
% -------------------------------------------------------------------------
cfo_estimation_mpower = 4;
cfo_estimation_sample_delay_block = 1; 
cfo_estimation_sample_delay_point = 2; 
cfo_estimation_block_length = 1000;
% Parameters of the phase differential algorithm (Leven)

% -------------------------------------------------------------------------
% Added CFO calculation
% -------------------------------------------------------------------------
tk = [0:1:length(symbs)-1]/symbol_rate;
% Sampling instants

phi0 = cfo_pp/2/cfo_var_freq;
pphi = phi0*cos(2*pi*cfo_var_freq*tk);
% Time varying phase, representing the sinusoidal variations of the CFO
% around its mean frequency

cfo_added = cfo_mean + cfo_var_freq*phi0*sin(2*pi*cfo_var_freq*tk);
% Calculation of the added CFO, for representation purpose

% -------------------------------------------------------------------------
% Add CFO
% -------------------------------------------------------------------------
symbs = symbs.*exp(1i*2*pi*cfo_mean*tk).*exp(-1i*pphi);
% CFO addition

% -------------------------------------------------------------------------
% Calculate delay of the estimator etc.
% Is actually also returned by updated version of function
% -------------------------------------------------------------------------
cfo_estimation_delay_block = (cfo_estimation_sample_delay_block + cfo_estimation_block_length - 1)/2;
cfo_estimation_delay_point = (cfo_estimation_sample_delay_point + 1 - 1)/2;
% Estimation of CFO estimation delay, in number of samples

cfo_compensable_max_block = symbol_rate/2/cfo_estimation_sample_delay_block/cfo_estimation_mpower;
cfo_compensable_max_point = symbol_rate/2/cfo_estimation_sample_delay_point/cfo_estimation_mpower;
% Maximum compensable CFO, in Hz


fprintf('\n\n\n%s\n','Test of CFO estimation by Leven method')
fprintf('%s\n\n','======================================')


fprintf('%s\n','Point estimation:')
fprintf('%s\t\t%i\n','Delay between samples:',cfo_estimation_sample_delay_block)
fprintf('%s\t%i\n','Size of averaging block:',1)
fprintf('%s\t\t\t%1.1f\n','Estimator delay:',cfo_estimation_delay_point)
fprintf(1,'%s\t%6.6f\t%s\n\n','Maximum compensable CFO: ',cfo_compensable_max_point/1.0e6, 'MHz');

fprintf('%s\n','Block estimation:')
fprintf('%s\t\t%i\n','Delay between samples:',cfo_estimation_sample_delay_block)
fprintf('%s\t%i\n','Size of averaging block:',cfo_estimation_block_length)
fprintf('%s\t\t\t%1.1f\n','Estimator delay:',cfo_estimation_delay_block)
fprintf(1,'%s\t%6.6f\t%s\n\n','Maximum compensable CFO: ',cfo_compensable_max_block/1.0e6, 'MHz');


fprintf(1,'%s\t%6.6f\t%s\n','Added mean CFO: ',cfo_mean/1.0e6, 'MHz');
fprintf(1,'%s\t\t%6.6f\t%s\n\n','Added max CFO: ',max(cfo_added)/1.0e6, 'MHz')

% -------------------------------------------------------------------------
% Check constellation after CFO addition
% -------------------------------------------------------------------------
plot_constellation(symbs,'plain','Constellation after CFO addition');

% -------------------------------------------------------------------------
% Add AWGN
% -------------------------------------------------------------------------
symbs = add_awgn(symbs,esn0_db);
% Add white Gaussian noise
    
% -------------------------------------------------------------------------
% Check constellation after AWGN addition
% -------------------------------------------------------------------------
plot_constellation(symbs,'plain','Constellation after AWGN addition');

% -------------------------------------------------------------------------
% Make a copy of the signal to process
% -------------------------------------------------------------------------
symbs_rx = symbs;
% Make a copy of signal.

% -------------------------------------------------------------------------
% CFO estimation
% -------------------------------------------------------------------------
cfo_estimate_block = cfo_estimation_leven(symbs_rx,cfo_estimation_mpower,cfo_estimation_sample_delay_block,cfo_estimation_block_length);
% CFO estimation with block averaging
cfo_estimate_point = cfo_estimation_leven(symbs_rx,cfo_estimation_mpower,cfo_estimation_sample_delay_point,1);
% 1-point CFO estimation

% -------------------------------------------------------------------------
% Compare applied CFO with estimates
% -------------------------------------------------------------------------
cfo_added_delayed_block = dsp_delay(cfo_added,cfo_estimation_delay_block);
cfo_added_delayed_point = dsp_delay(cfo_added,cfo_estimation_delay_point);
% Delay the added CFO by the value of the delay of the estimator

cfo_estimate_error_block = cfo_estimate_block*symbol_rate - cfo_added_delayed_block;
cfo_estimate_error_point = cfo_estimate_point*symbol_rate - cfo_added_delayed_point;
% Calculate CFO estimation error, in Hz

figure('Name','CFO estimate')
plot(cfo_estimate_block*symbol_rate/1.0e6,'g-');
hold on
plot(cfo_estimate_point*symbol_rate/1.0e6,'b--');
plot(cfo_added/1.0e6,'r--')
hline(cfo_mean/1.0e6,'k.-');
xlabel('sample');
ylabel('CFO (MHz)');
legend('CFO estimate (block)','CFO estimate (point)','added CFO')

figure('Name','CFO estimate : estimator delay check')
hold on
plot(cfo_estimate_block*symbol_rate/1.0e6,'go-')
hold on
plot(cfo_estimate_point*symbol_rate/1.0e6,'b*--')
plot(cfo_added/1.0e6, 'rs--')
xlabel('sample')
ylabel('CFO (MHz)')
legend('CFO estimate (block)','CFO estimate (point)','added CFO')
ylim((cfo_mean + [-0.1 0.1])/1.0e6)
hline(cfo_mean/1.0e6,'r-')
grid on

figure('Name','CFO estimate - point')
subplot(2,1,1)
plot(cfo_estimate_point*symbol_rate/1.0e6,'b-');
hold on
plot(cfo_added_delayed_point/1.0e6,'r--')
xlabel('sample');
ylabel('CFO (MHz)');
legend('CFO estimate (point)','added CFO (delayed)')
ylim([min(cfo_added)*0.95 max(cfo_added)*1.05]/1.0e6)
subplot(2,1,2)
plot(cfo_estimate_error_point/1.0e3,'k-');
xlabel('sample');
ylabel('CFO error (kHz)');
ylim([min(cfo_estimate_error_point(2*cfo_estimation_delay_point+1:end)) max(cfo_estimate_error_point(2*cfo_estimation_delay_point+1:end))]/1.0e3)

figure('Name','CFO estimate - block')
subplot(2,1,1)
plot(cfo_estimate_block*symbol_rate/1.0e6,'g-');
hold on
plot(cfo_added_delayed_block/1.0e6,'r--')
xlabel('sample');
ylabel('CFO (MHz)');
legend('CFO estimate (block)','added CFO (delayed)')
ylim([min(cfo_added)*0.95 max(cfo_added)*1.05]/1.0e6)
subplot(2,1,2)
plot(cfo_estimate_error_block/1.0e3,'k-');
xlabel('sample');
ylabel('CFO error (kHz)');
ylim([min(cfo_estimate_error_block(2*cfo_estimation_delay_block+1:end)) max(cfo_estimate_error_block(2*cfo_estimation_delay_block+1:end))]/1.0e3)


%%
% -------------------------------------------------------------------------
% CFO compensation
% -------------------------------------------------------------------------
symbs_comp_block = dsp_delay(symbs,cfo_estimation_delay_block).*exp(-1j*2*pi*cfo_estimate_block.*dsp_delay(tk,cfo_estimation_delay_block)*symbol_rate);
symbs_comp_block = symbs_comp_block(2*cfo_estimation_delay_block + 1:end);

symbs_comp_point = dsp_delay(symbs,cfo_estimation_delay_point).*exp(-1j*2*pi*cfo_estimate_point.*dsp_delay(tk,cfo_estimation_delay_point)*symbol_rate);
symbs_comp_point = symbs_comp_point(2*cfo_estimation_delay_point + 1:end);

symbs_original_delayed_point = dsp_delay(symbs_ref,cfo_estimation_delay_point);
symbs_original_delayed_point = symbs_original_delayed_point(2*cfo_estimation_delay_point + 1:end);

symbs_original_delayed_block = dsp_delay(symbs_ref,cfo_estimation_delay_block);
symbs_original_delayed_block = symbs_original_delayed_block(2*cfo_estimation_delay_block + 1:end);

% -------------------------------------------------------------------------
% Check constellation after CFO compensation
% -------------------------------------------------------------------------
plot_constellation(symbs_comp_block,'plain','Constellation after CFO compensation - block');
plot_constellation(symbs_comp_point,'plain','Constellation after CFO compensation - point');

% -------------------------------------------------------------------------
% Check symbols after CFO addition and compensation
% -------------------------------------------------------------------------
check_symbols_range = [1 5000];

figure('Name','Original and noisy symbols')
stem(angle(symbs_ref)/pi,'b');
hold on
stem(angle(symbs_rx)/pi,'rs')
xlabel('sample');
ylabel('phase (\times \pi rad)');
legend('original','after CFO + AWGN')
xlim(check_symbols_range)

figure('Name','Original and compensated symbols - point')
stem(angle(symbs_original_delayed_point)/pi,'b');
hold on
stem(angle(symbs_comp_point)/pi,'rs')
xlabel('sample');
ylabel('phase (\times \pi rad)');
legend('original','after CFO compensation: 1 point')
xlim(check_symbols_range)

figure('Name','Original and compensated symbols - block')
stem(angle(symbs_original_delayed_block)/pi,'b');
hold on
stem(angle(symbs_comp_block)/pi,'rs')
xlabel('sample');
ylabel('phase (\times \pi rad)');
legend('original','after CFO compensation: block')
xlim(check_symbols_range)

% -------------------------------------------------------------------------
% Figure alignment
% -------------------------------------------------------------------------
% Requires alignfigs function
% Available at https://github.com/nickhale/alignfigs
alignfigs(4)



