% -------------------------------------------------------------------------
% Test of diffenc_qam.m and diffdec_qam.m functions.
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
% QPSK
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Binary data stream
% -------------------------------------------------------------------------
nsymbols = 2^10;
% Number of symbols. Should be power of 2.

m = 4;
bits = generate_binary(nsymbols,m);
% Generate binary data

[words_dec,words_bin] = conv_bin2dec(bits,log2(m));
% Convert binary data into 2-bit words


% -------------------------------------------------------------------------
% Differential encoding to QPSK
% -------------------------------------------------------------------------
symbs_qpsk = diffenc_qpsk(words_dec);
% Differential encoding with the diffenc_qpsk function
symbs_4qam = diffenc_qam(bits,m);
% Differential encoding with the diffenc_qam function

plot_constellation(symbs_qpsk,'plain','QPSK - After differential encoding');
% Plot constellation after noise addition

fprintf('\n\n%s\n','Comparison between symbols generated by diffenc_qpsk and diffenc_qam');
symbs_comp = abs(symbs_4qam - symbs_qpsk);

mean_square_error_qpsk1_qam = sum(symbs_comp.^2)/nsymbols

% -------------------------------------------------------------------------
% WGN addition
% -------------------------------------------------------------------------
esn0_db = 10*4;
symbs_rx = add_awgn(symbs_qpsk,esn0_db);
% Add AWGN

plot_constellation(symbs_rx,'plain','QPSK - After noise addition');
% Plot constellation after noise addition




% -------------------------------------------------------------------------
% Differential decoding of QPSK
% -------------------------------------------------------------------------

[bits_qpsk_hard,symbs_diff_qpsk_hard] = diffdec_qpsk(symbs_rx,'hard');
[bits_qpsk_soft,symbs_diff_qpsk_soft] = diffdec_qpsk(symbs_rx,'soft');
[bits_4qam_hard,symbs_diff_4qam_hard] = diffdec_qam(symbs_rx,m);
% Differential decoding

plot_constellation(symbs_diff_qpsk_hard,'plain','QPSK - Differential constellation - QPSK hard decision');
plot_constellation(symbs_diff_qpsk_soft,'plain','QPSK - Differential constellation - QPSK soft decision');
plot_constellation(symbs_diff_4qam_hard,'plain','QPSK - Differential constellation - 4QAM hard decision');

% -------------------------------------------------------------------------
% BER calculation
% -------------------------------------------------------------------------
ber_qpsk_hard = calc_ber(bits_qpsk_hard,bits(3:end))
ber_qpsk_soft = calc_ber(bits_qpsk_soft,bits(3:end))
ber_4qam_hard = calc_ber(bits_4qam_hard,bits(3:end))
% BER calculation.



%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% 16-QAM
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Binary data stream
% -------------------------------------------------------------------------
nsymbols = 2^20;
% Number of symbols. Should be power of 2.

m = 16;
bits = generate_binary(nsymbols,m);
% Generate binary data

% -------------------------------------------------------------------------
% Differential encoding to 16-QAM
% -------------------------------------------------------------------------
symbs_16qam = diffenc_qam(bits,m);
% Differential encoding with the diffenc_qam function

plot_constellation(symbs_16qam,'plain','16-QAM - After differential encoding');
% Plot constellation after noise addition


% -------------------------------------------------------------------------
% WGN addition
% -------------------------------------------------------------------------
esn0_db = 40;
symbs_rx = add_awgn(symbs_16qam,esn0_db);
% Add AWGN.

plot_constellation(symbs_rx,'plain','16-QAM - After noise addition');
% Plot constellation after noise addition



% -------------------------------------------------------------------------
% Differential decoding of 16-QAM
% -------------------------------------------------------------------------
[bits_16qam_hard,symbs_diff_16qam_hard] = diffdec_qam(symbs_rx,m);
% Differential decoding


plot_constellation(symbs_diff_16qam_hard,'plain','16-QAM - Differential constellation - hard decision');




% -------------------------------------------------------------------------
% BER calculation
% -------------------------------------------------------------------------
ber_16qam_hard = calc_ber(bits_16qam_hard,bits(log2(m)+1:end))




%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% 64-QAM
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Binary data stream
% -------------------------------------------------------------------------
nsymbols = 2^10;
% Number of symbols. Should be power of 2.

m = 64;
bits = generate_binary(nsymbols,m);
% Generate binary data

% -------------------------------------------------------------------------
% Differential encoding to 64-QAM
% -------------------------------------------------------------------------
symbs_64qam = diffenc_qam(bits,m);
% Differential encoding with the diffenc_qam function

plot_constellation(symbs_64qam,'plain','64-QAM - After differential encoding');
% Plot constellation after noise addition


% -------------------------------------------------------------------------
% WGN addition
% -------------------------------------------------------------------------
esn0_db = 50;
symbs_rx = add_awgn(symbs_64qam,esn0_db);
% Add AWGN

plot_constellation(symbs_rx,'plain','64-QAM - After noise addition');
% Plot constellation after noise addition



% -------------------------------------------------------------------------
% Differential decoding of 64-QAM
% -------------------------------------------------------------------------
[bits_64qam_hard,symbs_diff_64qam_hard] = diffdec_qam(symbs_rx,m);
% Differential decoding. 


plot_constellation(symbs_diff_64qam_hard,'plain','64-QAM - Differential constellation - hard decision');


% -------------------------------------------------------------------------
% BER calculation
% -------------------------------------------------------------------------
ber_64qam_hard = calc_ber(bits_64qam_hard,bits(log2(m)+1:end))



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
























