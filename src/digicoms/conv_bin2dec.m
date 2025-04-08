function [words_dec,words_bin] = conv_bin2dec(bits_bin,nob)
% Converts vector of bits (binary) to vector of words (decimal)
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function converts a vector of bits to a vector of words of nob bits, 
% expressed in decimal or binary formats
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% [words_dec,words_bin] = conv_bin2dec(bits_bin,log2(m)); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% bits_bin          vector of bits [binary vector]
%
% nob               size of the words to generate, in number of bits
%                       [integer scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% words_dec         words, in decimal representation [integer vector]
%
% words_bin         binary representation of words [binary matrix]
%
%                       Matrix where the number of lines is the number of 
%                       words and the number of columns is nob.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

bits_bin = bits_bin(:)';
% Ensure the input bits are under the form of a line vector

words_bin = reshape(bits_bin,[nob length(bits_bin)/nob])';
% Matrix containing the bits of the binary representation of the symbols
% The number of lines corresponds to the number of symbols 
% The number of columns corresponds to the number of bits per symbol, i.e.
% log2(M).

words_dec = bin2dec(num2str(words_bin))';
% Decimal word vector

% words_bin_char = dec2bin(words_dec);
% Vector containing the binary representations of the symbols (char).


end