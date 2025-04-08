function symbs = mapping(words_dec,constellation)
% Mapping decimal words to complex symbols
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function does very little... It has just been created to be able to
% call a function to perform the mapping operation.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% symbs = mapping(words_dec,constellation);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% words_dec         decimal representations of the word stream to be
%                       encoded [integer vector]
%
% constellation     considered constellation [complex vector]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% symbs             encoded symbols [complex vector]
%
%                       The elements of symbs take values within the
%                       elements of constellation.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

symbs = constellation(words_dec + 1);
% That's it...

end