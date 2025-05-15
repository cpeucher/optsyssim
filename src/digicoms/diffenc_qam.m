function symbs = diffenc_qam(bits,m)
% Generation of differentially encoded symbols for square m-QAM modulation
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function differentially encodes a bit stream onto an m-QAM
% constellation. It divides the M-QAM constellation into a number of
% auxiliary QPSK subconstellations (1 for QPSK, 1 + 4 for 16-QAM, 1 + 4 +
% 16 for 64-QAM, etc) according to the principle described in e.g.
% J.-K. Hwang, et al., “Angle differential-QAM scheme for resolving phase 
% ambiguity in continuous transmission system,” International Journal of 
% Communication Systems, vol. 21, no. 6, pp. 631–641, Jun. 2008, 
% doi: 10.1002/dac.914.
%
% For instance a symbol of a 16-QAM constellation can be represented by
% S[k] = C[k] + D[k] where C[k] describes the quadrant centers, whose
% values are {2+2j, -2+2j, -2-2j, 2-2j} and D[k] describes the
% displacement with respect to the relevant center, whose values are {1+1j,
% -1+1j, -1-1j, 1-1j}.
% The quadrant is updated according to the first two bits in the 4-bit word
% representing the data to transmit in the 16-QAM constellation, and
% the displacement with respect to new quadrant center is updated according
% to the last two bits. 
% The scheme can be extended in a straightforward way to higher-order
% constellations (e.g. S[k] = B[k] + C[k] + D[k] for 64-QAM).
% Differential QPSK encoding is performed for each of the auxiliary
% constellations B, C, D etc.
% The differential encoding rule is as follows:
% ---------------------------------------------
%  word_dec  |   word_bin     |    phase shift
% ---------------------------------------------
%      0     |     00         |      0
%      1     |     01         |      pi/2
%      2     |     10         |      3*pi/2
%      3     |     11         |      pi
% ---------------------------------------------
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% symbs = diffenc_qam(bits,m);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% bits              bit stream to be encoded [integer 0/1 vector]
%
% m                 order of the modulation format [integer]
%
%                       To generate m-QAM constellation.
%                       Only square constellations are generated with this
%                       function, i.e. m = 4, 16, 64...
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% symbs             encoded symbols [complex vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

nbits = length(bits);
% Number of bits to encode

if rem(nbits,log2(m)) ~= 0
    % In case the number of bits provided does not correspond to an integer
    % number of symbols
    if nbits > log2(m)
        % If there is at least 1 symbol
        bits = bits(1:log2(m)*floor(nbits/log2(m)));
        % We truncate the bit stream so that it corresponds to an integer number of symbols
        fprintf('\n\n%s\n','diffenc_qam: the input bit stream has been truncated to an integer number of symbols.')
        % We warn the user.
    else
        error('diffenc_qam insufficient number of bits.')
        % We have fewer bits than necessary to encode one symbol.
    end
end
% Check input consistency. Normally this is not needed since we should
% generate a bit stream that is compatible with the chosen order of
% modulation.
    
[words_dec,~] = conv_bin2dec(bits,2);
% We divide the input bit stream in blocks of 2 bits (for encoding into
% QPSK sub-constellations), regardless of the order m of the modulation
% format.

nwords = length(words_dec);
% Number of 2-bit words

nwords_per_symbol = log2(m)/2;
% Number of 2-bit words necessary to encode one symbol

nsymbols = nwords/nwords_per_symbol;
% Number of symbols

nauxconst = log2(m)/2;
% Number of required auxiliary QPSK constellations

iquad = zeros(nauxconst,nsymbols);
% Initialise variable that will contains the indices of the quadrant of the
% auxilliary constellations
% Each line corresponds to an auxiliary constellation.
% Each column corresponds to a symbol, the first column being the
% arbitrarily set initial state, which will be added in the next few lines.
iquad0 = zeros(nauxconst,1);
% iquad(:,1) will contain the initial state = [0;0;...;0].
iquad = [iquad0,iquad];
% quad is a matrix of nauxconst lines and nsymbols + 1 columns. 
% It will later be truncated to nsymbols columns by removing the first
% column, which corresponds to the arbitrary initial state.
    

for isymbol = 1:nsymbols
    % Loop over symbols.
    for iauxconst = 1:nauxconst
        % Loop over auxiliary constellations.
        iquad(iauxconst,isymbol+1) = diffenc_qpsk_nextquadrant(iquad(iauxconst,isymbol),words_dec(1 + nwords_per_symbol * (isymbol - 1) + iauxconst - 1));
        % Determine the next quadrant index for each auxiliary
        % constellation.
    end
end
    
constellation_qpsk_std = [1+1i,-1+1i,-1-1i,1-1i];
% Standard QPSK constellation

symbs_aux = constellation_qpsk_std(iquad +1);
% Value of the auxilliary constellations calculated from the quadrant
% indices

const_mult = flipud(2.^([0:nauxconst-1]'));
% We need to multiply each subconstellation by the proper factor, i.e. 
% [1] for QPSK 
% [2;4] for 16-QAM
% [4;2:1] for 64-QAM
% etc.
% Here we prepare this column vector of multiplicative factors.  

symbs = sum(symbs_aux.*const_mult,1);
% We multiply the normalized subconstellations by the multiplicative factor
% and sum them symbol by symbol in order to obtain the final m-QAM symbol.

symbs = symbs(2:end);
% Remove the initial state.

end