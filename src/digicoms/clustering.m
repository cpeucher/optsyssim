function clust = clustering(symbs,words_dec,m)
% Data-aided clustering of received symbols
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function performs clustering of received symbols to a specified
% constellation (defined by the indices of the symbols in the
% constellation).
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% symbs_rx = normalise_constellation(symbs_rx,norm_es);
% cluster_da = clustering(symbs_rx,words_dec,m);
% % Data-aided clustering.
% symbs_cx = decision_qam_square_hard(symbs_rx,m);
% words_dec_rx = demapping(symbs_cx,constellation);
% cluster_dd = clustering(symbs_cx,words_dec_rx,m);
% % Decision-directed clustering.
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% symbs             complex symbols to cluster [complex vector]
%
% words_dec         decimal representation of the transmitted words
%                       [integer vector]
%
%                       words_dec can be a vector containing the binary
%                       representation of the original transmitted symbols,
%                       in which case we perform data-aided clustering, or
%                       words_dec can be the result of a decision, in
%                       which case we perform decision-directed clustering.
%
% m                 order of the modulation format [integer]
%
%                       Number of symbols in the constellation.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% clust             clusters of complex symbols [cell array]
%
%                       clust{ii,:} contains the elements (complex symbols)
%                       of the ii-th cluster, corresponding to transmitted
%                       word ii - 1 (decimal represention)
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

for ii = 1:m    
    clust{ii,:} = symbs(words_dec == ii-1);
end
% Separation into clusters based on the transmitted symbol (i.e. no
% thresholding is used here)
% clust{k,:} contains all the received samples corresponding to the 
% transmission of the symbol k -1 (in decimal form)

end