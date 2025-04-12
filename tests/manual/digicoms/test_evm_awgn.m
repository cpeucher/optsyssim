% -------------------------------------------------------------------------
% Test of error-vector magnitude (EVM) calculation over an additive white
% Gaussian noise (AWGN) channel.
% 
% Comparison between data-aided and decision directed clustering.
% For high signal-to-noise ratio (SNR), the results are identical. However, 
% when the SNR is reduced and the clusters start overlapping, the
% decision-directed and data-aided approaches start diverging.
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
% Global parameters
% -------------------------------------------------------------------------
nsymbols = 2^15;
% Number of symbols. Should be power of 2.

% -------------------------------------------------------------------------
% Define constellation
% -------------------------------------------------------------------------
m = 4;
% Constellation order. m points with log2(m) even (square constellation)
[constellation,norm_es,norm_emax] = define_constellation('qpsk_gray',m);
% Define QPSK constellation

% OR:

m = 16;
% Constellation order. M points with log2(M) even (square constellation)
[constellation,norm_es,norm_emax] = define_constellation('qam16_natural',m);
% Define 16-QAM constellation



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

plot_constellation(symbs,'plain','Original constellation');

scaling_factor = 100;
symbs_scaled = scaling_factor*symbs;
% To check that a scaling factor is irrelevant


% -------------------------------------------------------------------------
% AWGN
% -------------------------------------------------------------------------
esn0_db = 20;

% To see the difference between data-aided and decision-directed recovered
% constellation clusters, set
% esn0_db = 15;

symbs = add_awgn(symbs,esn0_db);
plot_constellation(symbs,'plain','After ASE addition');

symbs_scaled = add_awgn(symbs_scaled,esn0_db);
plot_constellation(symbs_scaled,'plain','After ASE addition - scaled',[-500:200:500]);


symbs_rx = symbs_scaled;
% Received symbols


% -------------------------------------------------------------------------
% Data-aided estimation of SNR
% -------------------------------------------------------------------------
for ii = 1:m    
    cluster_da{ii,:} = symbs_rx(words_dec == ii - 1);
end
% Separation into clusters based on the transmitteds symbol (i.e. no
% thresholding is used here).
% cluster_da{k,:} contains all the received samples corresponding to the 
% transmission of the symbol k (in decimal form). 


for ii = 1:m
    snr_cluster(ii) = mean(abs(cluster_da{ii,:}).^2)./(var(real(cluster_da{ii,:})) + var(imag(cluster_da{ii,:})));
end
% SNR per cluster

snr_retrieved_da_dB = 10*log10(mean(snr_cluster))
% Overall SNR
% This should correspond to the specified EsN0_dB.



% -------------------------------------------------------------------------
% Received signal normalisation: maximum energy
% -------------------------------------------------------------------------
for ii = 1:m       
    mean_constellation_da(ii) = mean(real(cluster_da{ii,:})) + 1i*mean(imag(cluster_da{ii,:}));
end
% We calculate the constellation of the mean points of each cluster (data
% aided) in the received signal.

norm_emax_rx = max(abs(mean_constellation_da).^2);
% Max energy over all mean symbols in the received constellation

norm_factor_emax = sqrt(norm_emax/norm_emax_rx);
% Normalisation factor towards a "standard constellation"

symbs_rx_norm_emax = symbs_rx * norm_factor_emax;
% The received constellation should now be rescaled to a standard
% constellation with clusters around {...,-3, -1, 1,3,...}

mean_constellation_da_norm_emax = mean_constellation_da * norm_factor_emax;
% We also normalise the mean constellation

% figure('Name','Normalised (Emax) mean retrieved constellation');
% scatter(real(meanConstellation_da_normEmax),imag(meanConstellation_da_normEmax));

plot_constellation(mean_constellation_da_norm_emax,'plain','Normalised (Emax - data-aided) mean retrieved constellation')


for ii =1:m
    cluster_da_norm{ii,:} = cluster_da{ii,:} * norm_factor_emax;
end

% figure('Name','Normalised (Emax) retrieved constellation');
% scatter(real(symbs_rx_normEmax),imag(symbs_rx_normEmax),'b','filled');
% hold on
% scatter(real(constellation),imag(constellation),'r','filled');
% maxI = max(abs(real(symbs_rx_normEmax)));
% maxQ = max(abs(imag(symbs_rx_normEmax)));
% maxIQ = max(maxI,maxQ);
% xlim(1.1*[-maxIQ maxIQ]);
% ylim(1.1*[-maxIQ maxIQ]);
% axis square
% xticks([-9:2:9]);
% yticks([-9:2:9]);


% figure('Name','Normalised (Emax) retrieved constellation - clusters');
% hold on;
% for ii =1:M
%     scatter(real(cluster_da_norm{ii,:}),imag(cluster_da_norm{ii,:}),'filled');
% end

% figure('Name','Normalised (Emax) retrieved constellation - heatmap');
% scatplot(real(symbs_rx_normEmax), imag(symbs_rx_normEmax)); 


plot_constellation(symbs_rx_norm_emax,'plain','Normalised (Emax - data-aided) retrieved constellation');
hold on;
scatter(real(constellation),imag(constellation),'r','filled');

plot_constellation(symbs_rx_norm_emax,'heat','Normalised (Emax - data-aided) retrieved constellation - heatmap');

plot_constellation(cluster_da_norm,'cluster','Normalised (Emax - data-aided) retrieved constellation - clusters');

%%

% -------------------------------------------------------------------------
% Received signal normalisation: energy per symbol
% -------------------------------------------------------------------------

norm_es_rx = mean(abs(symbs_rx).^2);
% Energy per symbol of the received (rescaled + noisy) constellation

norm_factor_es = sqrt(norm_es/norm_es_rx);
% Normalisation factor towards a "standard constellation"

symbs_rx_norm_es = symbs_rx * norm_factor_es;
% The received constellation should now be rescaled to a standard
% constellation with clusters around {...,-3, -1, 1,3,...}


% figure('Name','Normalised (Es) retrieved constellation');
% scatter(real(symbs_rx_normEs),imag(symbs_rx_normEs),'b','filled');
% hold on
% scatter(real(constellation),imag(constellation),'r','filled');

plot_constellation(symbs_rx_norm_es,'plain','Normalised (Es) retrieved constellation');
hold on;
scatter(real(constellation),imag(constellation),'r','filled');


symbs_cx = decision_qam_square_hard(symbs_rx_norm_es,m);
% Hard decision for square QAM
words_dec_rx = demapping(symbs_cx,constellation);
% Demapping


for ii = 1:m    
    cluster_dd_norm{ii,:} = symbs_rx_norm_es(words_dec_rx == ii - 1);
end
% Separation into clusters based on hard decision 
% cluster_dd_norm{k,:} contains all the received samples corresponding to the 
% transmission of the symbol k (in decimal form). 

% figure('Name','Normalised (Es) retrieved constellation - clusters');
% hold on;
% for ii =1:M
%     scatter(real(cluster_dd_norm{ii,:}),imag(cluster_dd_norm{ii,:}),'filled');
% end

plot_constellation(cluster_dd_norm,'cluster','Normalised (Es - decision-directed) retrieved constellation - clusters');


% -------------------------------------------------------------------------
% EVM calculation: data-aided
% -------------------------------------------------------------------------
for ii = 1:m          
    sum_ii(ii) = sum(abs(cluster_da_norm{ii,:} - constellation(ii)).^2);       
end

evm_max_da = sqrt(sum(sum_ii)/nsymbols)/max(abs(constellation));
evm_rms_da = sqrt(sum(sum_ii)/nsymbols)/ sqrt(mean(abs(constellation).^2));


% -------------------------------------------------------------------------
% EVM calculation: decision-directed
% -------------------------------------------------------------------------
for ii = 1:m        
    sum_ii(ii) = sum(abs(cluster_dd_norm{ii,:} - constellation(ii)).^2);     
end

evm_max_dd = sqrt(sum(sum_ii)/nsymbols)/max(abs(constellation));
evm_rms_dd = sqrt(sum(sum_ii)/nsymbols)/ sqrt(mean(abs(constellation).^2));



% -------------------------------------------------------------------------
% Display result
% -------------------------------------------------------------------------

fprintf('\n\n%s','EVM Test');
fprintf('\n%s','--------');

fprintf('\n\n%s\t%3.4f','EVM_rms (data-aided)=',evm_rms_da);
fprintf('\n%s\t%3.4f','EVM_max (data-aided) =',evm_max_da);
fprintf('\n\n%s\t%3.4f','EVM_rms (decision-directed)=',evm_rms_dd);
fprintf('\n%s\t%3.4f','EVM_max (decision-directed) =',evm_max_dd);


fprintf('\n\n%s\t%3.4f','Expected EVM_rms from Es/N0=',1/sqrt(10^(esn0_db/10)));

fprintf('\n\n%s\t%3.4f','k=EVM_rms/EVM_max (data-aided)=',evm_rms_da/evm_max_da);
fprintf('\n%s\t%3.4f','k=EVM_rms/EVM_max (decision-directed)=',evm_rms_da/evm_max_da);
fprintf('\n%s\t%3.4f','Expected value for PSK: k=',1);
fprintf('\n%s\t%3.4f','Expected value for QAM16: k=',sqrt(9/5));
fprintf('\n%s\t%3.4f','Expected value for QAM32: k=',sqrt(17/10));
fprintf('\n%s\t%3.4f\n\n','Expected value for QAM64: k=',sqrt(7/3));


[evm_max_dd_check,evm_rms_dd_check] = calc_evm(constellation,cluster_dd_norm);
[evm_max_da_check,evm_rms_da_check] = calc_evm(constellation,cluster_da_norm);




% -------------------------------------------------------------------------
% Figure alignment
% -------------------------------------------------------------------------
% Requires alignfigs function
% Available at https://github.com/nickhale/alignfigs
alignfigs(4)


% -------------------------------------------------------------------------
% Save data
% -------------------------------------------------------------------------
% file_name_data = [file_name_core_data '_' time_stamp]
% save(file_name_data,'');

