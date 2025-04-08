function [evm_max,evm_rms] = calc_evm(constellation,clust)
% Error vector magnitude (EVM) calculation
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the error vector magnitude of a received
% constellation (after clustring). Two normalisations (maximum energy and 
% energy-per-symbol of the constellation) are performed.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% [evm_max,evm_rms] = calc_evm(constellation,clust);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% constellation     standard constellation [complex vector]
%
% clust             clusters of complex symbols normalised to the
%                       standard constellation [cell array]
%
%                       clust{ii,:} contains the elements 
%                       (complex symbols) of the ii-th cluster 
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% evm_max           evm normalised by the maximum energy of the
%                       constellation [real scalar]
%
% evm_rms           evm normalised by the energy-per-symbol of the
%                       constellation [real scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

m = length(constellation);
% Number of elements in the constellation

nsymbols = 0;
% Initialise number of symbols

for ii = 1:m
    sum_ii(ii) = sum(abs(clust{ii,:} - constellation(ii)).^2);
    % Mean squared error
    nsymbols = nsymbols + length(clust{ii,:});
    % Update number of symbols
end

evm_rms = sqrt(sum(sum_ii)/nsymbols)/ sqrt(mean(abs(constellation).^2));
% EVM normalised by (average) energy per symbol of the constellation
evm_max = sqrt(sum(sum_ii)/nsymbols)/max(abs(constellation));
% EVM normalised by max constellation energy

end