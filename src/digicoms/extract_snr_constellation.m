function snr_db = extract_snr_constellation(clust)
% Extract the SNR from data-aided clusters of a constellation
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function extract the signal-to-noise ratio from data-aided clusters
% of a constellation.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% snr_db = extract_snr_constellation(cluster_da);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% clust             clusters of complex symbols [cell array]
%                       clust{ii,:} contains the elements (complex symbols)
%                       of the ii-th cluster, corresponding to transmitted
%                       word ii - 1 (decimal representation)
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% snr_db             SNR in dB [real scalar]
%                       This should correspond to Es/N0 (in dB).
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

for ii = 1:numel(clust)
    snr(ii) = mean(abs(clust{ii,:}).^2)./(var(real(clust{ii,:})) + var(imag(clust{ii,:})));
end
% SNR per cluster 

snr_db = 10*log10(mean(snr));
% Overall SNR
% This should correspond to the specified esn0_db

end