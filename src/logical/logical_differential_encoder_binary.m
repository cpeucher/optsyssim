function diffenc = logical_differential_encoder_binary(data)
% Differential encoder for DPSK or duobinary signal generation
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function performs differential encoding of a binary sequence. It is
% to be used for differential encoding in e.g. DPSK and optical duobinary
% transmitters.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% bit_pattern_diffenc = logical_differential_encoder_binary(bit_pattern);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% data              binary sequence to be differentially encoded
%                       [binary vector]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% diffenc           differentially encoded binary sequence
%                       [binary vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

nbits = length(data);
% Number of bits in the input sequence

diffenc = zeros(1,nbits);
% Vector where the differentially encoded sequence will be stored

diffenc(1) = xor(data(1),1);

for i = 2:nbits
    diffenc(i) = xor(diffenc(i-1),data(i));
end

end