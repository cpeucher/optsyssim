function bits = generate_binary(nsymbols,m)
% Generation of uniformly distributed binary random data
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates a stream of uniformly distributed binary random 
% data.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% bits_bin = generate_binary(nsymbols,m);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% nsymbols          number of symbols to generate [integer scalar]
%
%                       Very often a power of 2.
%
% m                 number of points in the constellation [integer scalar]
%
%                       The number of bits per symbol is log2(m).
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% bits              bit stream [binary vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

bits = logical(randi([0,1],1,nsymbols*log2(m)));

end