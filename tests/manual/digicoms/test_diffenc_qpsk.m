% -------------------------------------------------------------------------
% Test of dsp_qpsk_diffenc.m and dsp_qpsk_diffdec.m functions.
% 
% 1) We perform differential encoding through the 2 methods implemented in
% the diffenc_qpsk.m function. We highlight the round-off error 
% accumulation in 'method2'.
% We also compare with differential encoding performed in the special case
% of 4-QAM by the diffenc_qam.m function.
% We check that we can perform soft-decision compatible differential
% decoding in the absence of additive noise.
% 2) We illustrate the difference between the hard-decision differential
% decoding and the soft-decision compatible differential decoding. In the
% first approach, we start by performing a hard quadrant decision before
% determining the differential phases / corresponding 2-word bits. In the
% second approach, we calculate the differential phase / corresponding
% 2-word bits directly on the received samples, without earlier hard
% decision.
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

%--------------------------------------------------------------------------
% Switches
%--------------------------------------------------------------------------
do_debug = 1;
do_print = 0;


% -------------------------------------------------------------------------
% Figure settings
% -------------------------------------------------------------------------
fig.interpreter = 'latex';
fig.font_size = 18;

mred = [178 0 24]/255;

% line_color = linspecer(6,'qualitative');
% Use plot(x,y,'Color',line_color(1,:));



% -------------------------------------------------------------------------
% Global parameters
% -------------------------------------------------------------------------
nsymbols = 2^20;
% Number of symbols. Should be power of 2


% -------------------------------------------------------------------------
% Binary data stream
% -------------------------------------------------------------------------
fprintf('\n\n%s\n','Generating binary data');
tic
m = 4;
data = generate_binary(nsymbols,m);
% Generate binary data
toc

fprintf('\n\n%s\n','Converting binary data into m-bit words');
tic
[words_dec,words_bin] = conv_bin2dec(data,log2(m));
% Convert binary data into 2-bit words
toc

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Differential encoding: comparison of two methods.
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
fprintf('\n\n%s\n','Differential encoding: QPSK method1');
tic
symbs = diffenc_qpsk(words_dec,'method1');
toc

fprintf('\n\n%s\n','Differential encoding: QPSK method2');
tic
symbs_method2 = diffenc_qpsk(words_dec,'method2');
toc

fprintf('\n\n%s\n','Differential encoding: 4-QAM');
tic
symbs_4qam = diffenc_qam(data,m);
toc

fprintf('\n\n%s\n','Comparison between symbols generated by the 2 QPSK methods');
symbs_comp_qpsk = abs(symbs - symbs_method2);

mean_square_error_qpsk12 = sum(symbs_comp_qpsk.^2)/nsymbols

fprintf('\n\n%s\n','Comparison between symbols generated by the QPSK method1 and 4-QAM');
symbs_comp = abs(symbs_4qam - symbs);

mean_square_error_qpsk1_qam = sum(symbs_comp.^2)/nsymbols


figure('Name','error between the 2 QPSK methods')
plot(symbs_comp_qpsk)
xlabel('sample')
ylabel('|error|')

figure('Name','error between the QPSK method1 and 4-QAM')
plot(symbs_comp)
xlabel('sample')
ylabel('|error|')





% -------------------------------------------------------------------------
% Check constellation and BER
% -------------------------------------------------------------------------
plot_constellation(symbs,'plain','Differentially-encoded constellation');

[data_rx_soft,~] = diffdec_qpsk(symbs,'soft');
[data_rx_hard,~] = diffdec_qpsk(symbs,'hard');
% Differential decoding. 

ber_soft = calc_ber(data_rx_soft,data(3:end))
ber_hard = calc_ber(data_rx_hard,data(3:end))
% BER calculation.



%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Illustration of the difference between the two differential decoding
% schemes.
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Add white Gaussian noise
% -------------------------------------------------------------------------
esn0_db = 20;
symbs_rx = add_awgn(symbs,esn0_db);
% Add AWGN

plot_constellation(symbs_rx,'plain','After noise addition');
% Plot constellation after noise addition


% -------------------------------------------------------------------------
% Differential decoding: hard decision
% -------------------------------------------------------------------------
[bits_hd,symbs_diff_hd] = diffdec_qpsk(symbs_rx,'hard');
% Differential decoding with hard decision

plot_constellation(symbs_diff_hd,'plain','differential constellation after HD');
% Plot differential constellation after hard decision

ber_hd = calc_ber(bits_hd,data(3:end))
% BER calculation for soft decision

% -------------------------------------------------------------------------
% Differential decoding: soft decision
% -------------------------------------------------------------------------
[bits_sd,symbs_diff_sd] = diffdec_qpsk(symbs_rx,'soft');
% Differential decoding with soft decision

plot_constellation(symbs_diff_sd,'plain','differential constellation after HD');
% Plot differential constellation after hard decision

ber_sd = calc_ber(bits_sd,data(3:end))
% BER calculation for soft decision



% -------------------------------------------------------------------------
% Figure alignment
% -------------------------------------------------------------------------
% Requires alignfigs function
% Available at https://github.com/nickhale/alignfigs
% alignfigs(2)


% -------------------------------------------------------------------------
% Save data
% -------------------------------------------------------------------------
% file_name_data = [file_name_core_data '_' time_stamp]
% save(file_name_data,'');

