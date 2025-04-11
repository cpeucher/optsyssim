function symbs = diffenc_qpsk(words_dec,method)
% Generation of differentially encoded symbols for QPSK modulation
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates normalized differentially encoded QPSK symbols
% from a binary word stream, where the words take values within {0,1,2,3}.
% Gray encoding is used, meaning that:
% if words_dec(i) = 0 (binary '00'), then phase(i) = phase(i-1)
% if words_dec(i) = 1 (binary '01'), then phase(i) = phase(i-1) + pi/2
% if words_dec(i) = 2 (binary '10'), then phase(i) = phase(i-1) + 3*pi/2
% if words_dec(i) = 3 (binary '11'), then phase(i) = phase(i-1) + pi
% The symbols are generated on the constellation [1+1i,-1+1i,-1-1i,1-1i].
% For the fun of it, 2 implementation methods are provided. 
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% symbs = diffenc_qpsk(words_dec);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% words_dec         decimal representations of the word stream to be 
%                       encoded [integer vector]
%
% method            implementation method [string]
%
%                       method = 'method1'
%                           Loop over the input words and definition of the
%                           quadrant index in the constellation.
%
%                       method = 'method2' 
%                           Do not use!
%                           For pedagogical purpose only.
%                           Search for the different values of the words
%                           and define the corresponding differential
%                           phase, then use cumsum to calculate the
%                           absolute phase.
%
%                       This is an optional parameter.
%                       If no input string is provided, then 'method1',
%                       which is somehow faster, is chosen.%
%                       Furthermore, round off errors are present in
%                       'method2' due to the use of cumsum of the phases.
%                       Consequently 'method2' is kept here just to
%                       illustrate that different approaches are possible,
%                       but the use of 'method1' should be favoured.
%                       This corresponds to the default behaviour of the
%                       function when no explicit method string is
%                       provided.
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

if nargin == 1    
    method = 'method1';    
end
% In case no method is explicitely stated, we choose 'method1', which
% appears to be faster and is not prone to round off errors.


switch method
    % Switch over calculation method
    
    case 'method1'
        
        constellation = [1+1i, -1+1i, -1-1i, 1-1i];
        % Define constellation
        % The index iquad in constellation(iquad+1) defines the quadrant:
        % iquad = [0,1,2,3]
        
        iquad = zeros(1,length(words_dec) + 1);
        % Pre-initialise quadrant index vector
        
        iquad(1) = 0;
        % We arbitrarily take the initial (i.e. previous) symbol in the 1st
        % quadrant
        
        for iword = 2:length(words_dec) + 1
            % Loop over words / qpsk sysmbols and assign quadrant to the next
            % symbol
            
            iquad(iword) = diffenc_qpsk_nextquadrant(iquad(iword - 1),words_dec(iword - 1));
            % Determine quadrant index based on 2-bit word value and
            % previous quadrant index
            
        end
        % End of loop over words / symbols
        % iquad is a vector of length nsymbols + 1, since we initialised it 
        % with the previous quadrant value.
        iquad = iquad(2:end);
        % We remove the initial value from the quadrant indices vector.
        
        symbs = constellation(iquad + 1);
        % Mapping
        
        
    case 'method2'        
        
        diff_phase = zeros(1,length(words_dec));
        % Preinitialise a differential phase vector        
        
        diff_phase(words_dec == 0) = 0;
        diff_phase(words_dec == 1) = pi/2;
        diff_phase(words_dec == 2) = 3*pi/2;
        diff_phase(words_dec == 3) = pi;
        % Assign the value of differential phase depending on the word to
        % be transmitted       
        
        symbs_phase = cumsum(diff_phase) + pi/4;
        % Calculate absolute phase from the differential phase. We add pi/4
        % since we assume, as for method1, that the previous symbol was
        % 1+1i (without loss of generality...)        
        
        symbs = sqrt(2)*exp(1i*symbs_phase);
        % Complex symbol (standard constellation) based on the provided
        % absolute phase value.
        
end
% End of switch over differential encoding method


end