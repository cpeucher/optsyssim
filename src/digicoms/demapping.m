function words_dec = demapping(symbs_cx,constellation)
% Demapping of digital symbols
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements demapping of digital symbols after decision to
% a standard constellation.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% words_dec = demapping(symbs_cx,constellation);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% symbs_cx          complex symbols after decision [complex vector]
%                       The elements of symbs_cx take value within the
%                       ensemble defining the constellation.
%
% constellation     definition of the considered constellation 
%                       [complex vector]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% words_dec         decimal representation of the word corresponding to
%                       the received symbols symbs_cx in the the
%                       constellation constellation [integer vector]
%
%                       This is index - 1 where
%                       constellation(index) = symbs_cx
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

words_dec = 0;
% Initialise

for isymb = 1:length(constellation)
    
    idx(isymb,:) = symbs_cx == constellation(isymb);
    
    words_dec = words_dec + idx(isymb,:) * isymb;
    
end

words_dec = words_dec - 1;
% In Matlab vectors are indexed from 1 to M, but for M-ary modulation, the
% decimal representation of the words are in 0 to M - 1

end