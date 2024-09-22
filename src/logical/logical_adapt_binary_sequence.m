function seq = logical_adapt_binary_sequence(seq,nsymbols)
% Adapt binary sequence to simulation time window
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function adapts the number of bits in a binary sequence to the the
% duration of the time window. It either truncates a too long pattern, or
% repeat (and possibly truncate) a too short pattern.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% seq = logical_adapt_binary_sequence(seq,nsymbols); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% seq               input binary sequence [binary vector]
%
% nsymbols          number of bits to which the sequence will be adapted
%                       [integer]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% seq               output binary sequence [binary vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% REMARKS:
% -------------------------------------------------------------------------
% Should also work for more general complex sequences. The output type is
% dictated by the input type.
%
% -------------------------------------------------------------------------
% TO DO:
% -------------------------------------------------------------------------
% Very old function. Pre Matlab literacy. Neither efficient nor elegant. 
% To clean up some day... 
% But has been doing the job so far...
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

pattern_length = length(seq);
% Length of the input sequence

k = log2(nsymbols/pattern_length);

if pattern_length > nsymbols
    seq = seq(1:nsymbols);
    % In case the pattern length is longer than the number of symbols
    % corresponding to the simulation time window, we truncate the pattern
    % and only keep the first Nsymbols bits.
    message = 'Binary pattern length is longer than time window duration: pattern will be truncated.';
    disp('logical_adapt_binary_sequence:');
    disp(message);
    
elseif pattern_length < nsymbols
    
    if mod(k,1) == 0
        % First we check if the total number of bits is equal to the
        % pattern length multiplied by a power of 2.
        for i=1:k
            seq = cat(2,seq,seq);
        end
    else
        % If the number of symbols is not a power of 2 multiple of the
        % pattern length:
        seq = seq(ones(1,ceil(nsymbols/pattern_length)),:)';
        seq = seq(:)';
        % The pattern is repeated ceil(nsymbols/pattern_length) times.
        seq = seq(1:nsymbols);
        % The repeated pattern is truncated so that its length is nsymbols.
    end
    message = 'Binary pattern length is shorter than time window duration: pattern will be repeated.';
    disp('logical_adapt_binary_sequence:');
    disp(message);
    
else
    % In case the pattern length matches the duration of the time window,
    % nothing needs to be done.
    % message = 'Binary pattern length matches time window duration.';
end

end