function symbs = normalise_constellation(symbs,norm_es)
% Normalisation of received constellation using energy per symbol
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function performs the normalisation of a received constellation 
% using its energy per symbol.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% symbs_rx = normalise_constellation(symbs_rx,norm_es);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% symbs             received symbols to normalise [complex vector]
%
% norm_es           energy per symbol of the normalised constellation
%                   corresponding to the received symbols [real scalar]
%
%                       This value is returned when generating the standard
%                       constellation with:
%                       [constellation,norm_es,norm_emax] = ...
%                           define_constellation(type,m);
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% symbs             symbols normalised to the standard constellation
%                       [complex vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

normEs_rx = mean(abs(symbs).^2);
% Energy per symbol of the received constellation

symbs = symbs * sqrt(norm_es/normEs_rx);
% The received constellation should now be rescaled to a standard
% constellation with clusters around {...,-3, -1, 1,3,...}

end