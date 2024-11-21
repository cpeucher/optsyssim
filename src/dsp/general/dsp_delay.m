function seqout = dsp_delay(seqin,sample_delay,initial_state)
% Digital delay
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements a digital delay on an input sequence. 
% If the sequence is input as a line vector, the delay is applied according
% to:
% sequin = [a1 a2 a3 a4 a5 a6 a7 a8 a9 a10]
% sample_delay = 3;
% initial_state = [i1 i2 i3];
% seqout = [i1 i2 i3 a1 a2 a3 a4 a5 a6 a7];
% The lengths of the input and output sequences are identical.
% The behavior of the function is compatible with that of dsp_fir_linear
% seqout = dsp_fir_linear(sequin,[0 0 0 1],[i3 i2 i1 i0])
% Example: 3 sample delay
% z = [0 -1 -2 -3];
% b = [0 0 0 1];
% y = dsp_fir_linear(x,b,z);
% and
% y = dsp_delay(x,3,[-2 -1 0]);
% return the same sequence.
% The function can act on several input sequences in parallel if seqin in
% a matrix where each line is a sequence.
% If no initial_state is provided, it will be set to a sequence of 0.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% a = linspace(1,10,10);
% b = dsp_delay(a,5);
% c = dsp_delay(a,5,[1 2 3 4 5]);
% aa = linspace(1,50,50);
% aa = reshape(aa,5,10);
% bb = dsp_delay(aa,3);
% cc = dsp_delay(aa,3,[9 10 11]);
% dd = dsp_delay(aa,3,[9 10 11;12 13 14;15 16 17;18 19 20;21 22 334]);
% 
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% seqin             input sequence [matrix]
%
%                       Each line of the matrix corresponds to an
%                       individual sequence to delay by the same amount of
%                       samples (parallel processing).
%                       A line vector will represent a single sequence.
%
% sample_delay      sample delay to apply to all signals in seqin [integer]
%
%                       sample_delay is a positive integer, i.e. we can
%                       only delay the input sequence, not advance it.
%
% initial_state     initial state of the buffer [matrix]
%
%                       If no initial_state is provided, it will by
%                       default be taken equal to [0 0 0] (for all
%                       sequences if seqin is a matrix).
%
%                       If seqin is a matrix and initial_state is a line
%                       vector, then the same initial_state will be applied
%                       to all the sequences (lines) in seqin.
%                       The length of initial_state should be equal to 
%                       sample_delay.
%
%                       If seqin is a matrix and initial_state is a matrix,
%                       then individual initial states will be applied to
%                       each sequence (line) in seqin. 
%                       The number of lines of initial_state should
%                       then be equal to the number of lines of seqin.
%                       The number of columns of initial_state should be
%                       equal to sample_delay.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% seqout             output sequence [matrix]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


[nl_seqin,nc_seqin] = size(seqin);
% Number of lines and columns in seqin. Will be used to determine whether a
% single sequence is processed, or parallel processing is performed on
% multiple sequences of the same length.    

switch nargin
    
    case 2
        % No initial state is provided. It will be taken equal to zero by
        % default.        
        initial_state = zeros(nl_seqin,sample_delay);
        
    case 3
        % An initial state is provided.        
        [nl_init,nc_init] = size(initial_state);
        
        if nc_init ~= sample_delay
            error('dsp_delay: the length of the initial sequence should be equal to the sample delay.');
        end
        
        if nl_init == 1
            % If the initial state is a line vector, the same values will
            % be applied to all sequences.
            initial_state = repmat(initial_state,[nl_seqin,1]);
        elseif nl_init ~= nl_seqin
            error('dsp_delay: the number of initial states should be equal to 1 or to the number of sequences.');
        end
        
    otherwise
        
        error('dsp_delay: too many input arguments.');
        
end
            
            
seqout= [initial_state,seqin];
seqout = seqout(:,1:nc_seqin);

end