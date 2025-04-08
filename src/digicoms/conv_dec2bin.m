function bits_bin = conv_dec2bin(words_dec,nob)
% Converts vector of words (decimal) to vector of bits (binary)
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function converts a vector of words (decimal) to a vector of bits
% (binary).
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% bits_bin = conv_dec2bin(words_dec,nob);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% words_dec         vector of numbers in decimal representation 
%                       [real vector]
%                       The numbers will be converted to their binary 
%                       representation with nob bits.
%
% nob               number of bits used for binary representation 
%                       [integer scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% bits_bin          vector of bits corresponding to the binary
%                       representation of the decimal numbers in 
%                       words_dec [binary vector]
%                       The length of the vector is equal to 
%                       length(words_dec)*nob.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

a = dec2bin(words_dec(:)',nob) - '0';

bits_bin = reshape(a',[1,numel(a)]);

bits_bin = logical(bits_bin);
% Convert the array to a logical array.


end