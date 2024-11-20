function [seq_i,seq_q,r,s,periodicity_check] = logical_differential_encoder_dqpsk(type,u,v,varargin)
% Differential encoder for DQPSK 
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function performs differential encoding of binary sequences for 
% use in optical DQPSK transmitters.
% The function can also be used to calculate the expected patterns at the 
% outputs of the differential receiver when arbitrary patterns are used to
% drive the modulators in the QPSK transmitter.
% The notations and equations are consistent with
% C. Peucheret, "Note on Optical DQPSK Modulation", Technical University
% of Denmark, September 2003 and, not really surprisingly, with 
% T. Tokle, "Optimised dispersion management and modulation formats for
% high speed optical communication systems", Ph.D. thesis, Technical
% University of Denmark, Dec. 2004.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% [bit_pattern_i,bit_pattern_q,r,s,periodicity_check] = logical_differential_encoder_dqpsk('parallel',bit_pattern1,bit_pattern2);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% type              transmitter type for which the differential encoded
%                       data will be generated [string]
%
%                       This relates to the generation of proper 
%                       differentially encoded sequences for 
%                       interferometric detection of optical DQPSK signals.
%
%                       type= 'serial': for 2 phase modulators (pi & pi/2), 
%                           or equivalently a 4-level electrical driving 
%                           signal with proper levels used to drive a
%                           single phase modulator
%
%                       type= 'parallel': for IQ modulator (dual parallel 
%                           Mach-Zehnder)
%
%                       type= 'serial_none': enables to calculate the 
%                           expected sequences after the differential '
%                           receiver when the u, and v data sequences are
%                           directly applied to the modulators at the
%                           serial transmitter, i.e. no differential 
%                           encoding of u and v is performed.
%
%                       type= 'parallel_none': enables to calculate the 
%                           expected sequences after the differential 
%                           receiver when the u, and v data sequences are
%                           directly applied to the modulators at the
%                           parallel transmitter, i.e. no differential 
%                           encoding of u and v is performed.
%
% u                 original binary data sequence [binary vector]
%
% v                 original binary data sequence [binary vector]
%
% varargin          varargin = [seq_i(1) seq_q(1)]
%                       [2-element binary vector]
%
%                       optional initial state of the differential encoder
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% seq_i             differentially encoded binary data sequence 
%                       [binary vector]
%
%                       Corresponds to in-phase in case a dual parallel 
%                       IQ modulator is used
%
% seq_q             differentially encoded binary data sequence 
%                       [binary vector]  
%
%                        Corresponds to quadrature in case a dual 
%                        parallel IQ modulator is used.
%
% r                  expected received pattern in the in-phase 
%                        differential receiver [binary vector]
%
% s                   expected received pattern in the quadrature 
%                       differential receiver [binary vector]
%
% periodicity_check   binary array specifying whether the I or Q 
%                       patterns are periodic with a period equal to
%                       their length [binary array]
%
%                           if periodicity_check(1) = 0, the I pattern is
%                               periodic.
%                           if periodicity_check(2) = 0, the Q pattern is 
%                               periodic.
%                           If periodicity_check = 1, make sure to ignore
%                           the first bit when comparing the received data
%                           to the original sequence. If the comparison is 
%                           made to r and s (for "standard" differential 
%                           receiver, then no problem is expected, since r 
%                           and s may differ from u and v by their first 
%                           bits.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% CREDITS:
% -------------------------------------------------------------------------
% Contains code by Torger Tokle
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

if length(u) ~= length(v)
    error('logical_differential_encoder_dqpsk: input data streams should have the same length.');
else
    nbits = length(u);
    % Number of bits in the input sequence
end
% Check that the input bit sequences u and v have the same lengths


optargin = size(varargin,2);
% Number of optional arguments

if optargin == 1
    init = varargin{1};  
else
    init = [0 1];
    % Arbitrarily initialise the differential encoder to the state [0 1].
    % Did not look at the theory, but does not seem to matter at all. If the
    % differentially encoded patterns are not Nbits periodic, then changing
    % the initial state will not change the status and will result in
    % either the first bit of r, the first bit of s, or both, to be wrong,
    % regardless of the choice of the initial state. So, there is no gain
    % in trying other initial states.
end
% Read optional initialisation state of the differential encoder   

type_tx = regexp(type,'\_','split');
type_tx = type_tx(1);
% Extract the type of transmitter from the type input field
% i.e.:
% 'Serial'   & 'Serial_None'   -> 'Serial'
% 'Parallel' & 'Parallel_None' -> 'Parallel'
        
u = [u u(1)];
v = [v v(1)];
% Extend the original data by adding the first bit to the end. u and v are
% now nbits+1 long.

[seq_i,seq_q] = logical_differential_encoder_dqpsk_calc_patterns(type,u,v,init(1),init(2));
% Calculate the differentially encoded sequences for nbits + 1.

periodicity_check = [abs(seq_i(nbits+1)-seq_i(1)) abs(seq_q(nbits+1)-seq_q(1))];
% Check the periodicity of the generated pattern by comparing the last and
% the first bits in the differentially encoded sequences

% global keep_log
% if keep_log.verbose == 1    
%     fprintf(1,'\n%s\n','Logical_DifferentialEncoder_DQPSK: ');
%     fprintf(1,'%s\t%i\n','I pattern periodicity error: ',periodicity_check(1));
%     fprintf(1,'%s\t%i\n','Q pattern periodicity error: ',periodicity_check(2));
%     fprintf(1,'%s\n\n','If periodicity error = 1, make sure to ignore the first bit when comparing the received data to the original sequence.');    
% end

seq_i = seq_i(1:nbits);
seq_q = seq_q(1:nbits);
% Restrict the differentially encoded sequences to their first nbits again

[r,s] = logical_differential_encoder_dqpsk_received_patterns(type_tx,seq_i,seq_q);
% Calculate the expected received patterns 


end
% End of main function logical_differential_encoder_dqpsk


% -------------------------------------------------------------------------
% Calculate differentially encoded patterns
% -------------------------------------------------------------------------
function [seq_i,seq_q] = logical_differential_encoder_dqpsk_calc_patterns(type,u,v,seq_i1,seq_q1)
% Differentially encoded DQPSK patterns for a given initial state of the
% encoder

seq_i = zeros(1,length(u));
seq_q = zeros(1,length(u));
% Arrays where the differentially encoded sequences will be stored

seq_i(1) = seq_i1;
seq_q(1) = seq_q1;
% Initialise the differential encoder

if strcmp(type,'parallel')
    
    for k= 2:length(u)
        seq_i(k) = not(xor(u(k),v(k))) & xor(u(k),seq_i(k-1)) | xor(u(k),v(k)) & xor(v(k),seq_q(k-1));
        seq_q(k) = not(xor(u(k),v(k))) & xor(v(k),seq_q(k-1)) | xor(u(k),v(k)) & xor(u(k),seq_i(k-1));
    end
    % Differential encoder for parallel transmitter
    
elseif strcmp(type,'serial')
    
    for k=2:length(u)
        seq_i(k) = not(xor(v(k),seq_i(k-1))) & not(seq_q(k-1)) | not(xor(u(k),seq_i(k-1))) & seq_q(k-1);
        seq_q(k) = xor(xor(u(k),v(k)),seq_q(k-1));
    end
    % Differential encoder for serial transmitter
    
elseif strcmp(type,'serial_none')
    
    seq_i = u;
    seq_q = v;
    % No differential encoder for serial transmitter; just used to
    % calculate the expected received pattern
    
elseif strmcp(type,'parallel_none')
    
    seq_i = u;
    seq_q = v;
    % No differential encoder for parallel transmitter; just used to
    % calculate the expected received pattern
    
else
    
    error('logical_differential_encoder_dqpsk_calc_patterns: differential encoder type not implemented.');
    
end

end
% End of logical_differential_encoder_dqpsk_calc_patterns function



% -------------------------------------------------------------------------
% Calculate received patterns after interferometric detection
% -------------------------------------------------------------------------
function [r,s] = logical_differential_encoder_dqpsk_received_patterns(type,seq_i,seq_q)
% Received patterns for DQPSK transmitter and interferometric detection

% Re-used code (c) Torger Tokle - 2006

if strcmp(type,'parallel')
    
    Phase = pi/4*(1 + 4*(seq_i + 0.5*xor(seq_i,seq_q)));
    % Phase generated by the dual-parallel Mach-Zehnder modulator
    
    DiffPhase = Phase-circshift(Phase',1)';
    % Differential phase
    
    r = logical(zeros(1,length(Phase)));
    s = logical(zeros(1,length(Phase)));
    % Initialise received arrays
    
    r = logical(round(0.5+0.5*(sin(DiffPhase)-cos(DiffPhase))));
    s = logical(round(0.5+0.5*(-sin(DiffPhase)-cos(DiffPhase))));
    % Received data after each balanced detector
    
elseif strcmp(type,'serial')
    
    Phase = pi*(seq_i + 0.5*seq_q);
    
    DiffPhase = Phase - circshift(Phase',1)';
    % Differential phase
    
    r = logical(zeros(1,length(Phase)));
    s = logical(zeros(1,length(Phase)));
    % Initialise received arrays
    
    r = logical(round(0.5+0.5*(sin(DiffPhase)-cos(DiffPhase))));
    s = logical(round(0.5+0.5*(-sin(DiffPhase)-cos(DiffPhase))));
    % Received data after each balanced detector
    
    r = ~r;
    s = ~s;
    % Logically invert data 
    
else
    
    error('logical_differential_encoder_dqpsk_received_patterns: DQPSK transmitter type not implemented.');
    
end

end
% End of logical_differential_encoder_dqpsk_received_patterns function

