function seq = logical_prbs(params)
% Generation of pseudo-random binary sequences
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function returns pseudo-random binary sequences obtained using shift
% registers.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_prbs.type = 'shift_register';%'de_bruijn';
% params_prbs.order = 7;
% params_prbs.poly = [7 6 0];%[7 3 0];[7 1 0];
% params_prbs.seed = [1 1 1 0 1 1 0];
% bit_pattern = logical_prbs(params_prbs);
% bit_pattern = logical_adapt_binary_sequence(bit_pattern,nsymbols);
%
% Other generating polynomials:
% params_prbs.order = 10;
% params_prbs.poly = [10 3 0];%[10 7 0];
% params_prbs.seed = [0 1 1 0 1 1 1 0 0 1];
% 
% params_prbs.order = 11;
% params_prbs.poly = [11 9 0];%[11 8 5 2 0];
% params_prbs.seed = [0 1 1 0 1 1 1 0 0 1 0];
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% params            PRBS generation parameters [structure]
%
%                       params.type
%                           sequence type [string]
%
%                           params.type = 'shift_register'
%                               standard shift register sequence of period 
%                               2^params.order -1
%                               An extra bit is added so that the length of 
%                               the returned sequence is 2^params.order. 
%                               This bit is set to the value of the first 
%                               bit of the sequence.
%
%                           params.type = 'de_bruijn' 
%                               2^params.order - 1 shift register sequence 
%                               to which a zero has been added at the end 
%                               of the run of consecutive zeros of 
%                               length params.order - 1
%                               Thus the 'de_bruijn' contains a
%                               run of 'params.order' consecutive zeros,
%                               unlike the 'shift_register'.
%
%                       params.order
%                           length of the shift register [integer]
%
%                           Use N to generate a 2^N length sequence
%
%                       params.poly
%                           generating polynomial [integer vector]
%
%                           The array contains the indices of the taps 
%                           use to generate the linear shift register 
%                           sequence.
%                           Use [7 6 0] for x^7+x^6+1=0, meaning that the 
%                           outputs of taps 7 and 6 are xor'ed and fed back 
%                           to the input of the shift register
%
%                       params.seed
%                           initial state of the shift register
%                           [binary vector]
%
%                           Should be a vector of length equal to 
%                           params.order. If seq=[a1 a2 a3 ... aN], where N 
%                           is the order of the PRBS, then a1 correponds to 
%                           the first bit of the generated sequence.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% seq               generated shift-register sequence [binary vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% REMARKS:
% -------------------------------------------------------------------------
% Pre-Matlab litteracy version. The PRBS could be generated without the
% loops.
%
% -------------------------------------------------------------------------
% TO DO:
% -------------------------------------------------------------------------
% TODO: implement faster algorithms to generate de Bruijn and shift
%       register sequences
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

sequence_length = 2^params.order;
% Length of the sequence to generate

prbs = zeros(1,sequence_length - 1);
% To store the shift register sequence of period 2^params.order - 1

seq = zeros(1,sequence_length,'logical');
% To store the returned sequence of length 2^params.order

% First we convert the generating polynomial to an array genpoly
% params.poly = [N m 0] means x^N + x^m + 1 = 0 is the generator 
% polynomial, which can then be converted to the Matlab format:
% genpoly = [1 0 0 1 0 ... 0 1]
%            |     |         |
%          deg N  deg m     deg 0
polyn = sort(params.poly,'descend');
% We ensure the taps numbers in the array poly are in decreasing order.
if polyn(1) ~= params.order
    % We check that the degree of the generating polynomial is consistent
    % with the order of the sequence.
    error('logical_prbs: PRBS order and degree of generating polynomial do not match.');
else
    genpoly = zeros(1,params.order + 1);
    % Initialise the generating polynomial to zero
    for i = 1:length(polyn)
        genpoly(params.order - polyn(i) + 1) = 1;
    end
end

if length(params.seed) ~= params.order
    % We check that the seed array has a length equal to order.
    error('logical_prbs: seed size does not match the requested PRBS order.');
else
    register = fliplr(params.seed);
    % Initialise the shift register
    for i = 1:sequence_length - 1
        prbs(i) = register(params.order);
        % Exit the last value of the present state of the shift register
        c = register(params.order)*genpoly(1);
        % Initialise the value of c that will be xor'ed and pushed back to
        % the input of the register. Hopefully genpoly(1) is equal to 1!
        for j = 1:params.order - 1
            c = xor(c,register(params.order - j)*genpoly(j + 1));
        end
        % Calculate next value to enter the register
        register = [c register(1:params.order - 1)];
        % Shift the content of the register to the right and insert the
        % calculated value
    end
end

% -------------------------------------------------------------------------
% Process the output sequence, depending on its type.
% -------------------------------------------------------------------------
if strcmp(params.type ,'de_bruijn')
    
    for i = 1:sequence_length - 1
        seq_run = zeros(1,params.order - 1);
        
        for k = i:i + params.order - 2
            if k > sequence_length - 1
                seq_run(k) = prbs(k - sequence_length + 1);
            else
                seq_run(k) = prbs(k);
            end
        end
        
        if sum(seq_run) == 0
            index = i + params.order - 2;
            % Return the index of the last zero of the run of order-1 zeros
        end
    end
    
    for i = 1:index
        seq(i) = prbs(i);
    end
    seq(index + 1) = 0;
    for i = index + 2:sequence_length
        seq(i) = prbs(i - 1);
    end
    
elseif strcmp(params.type ,'shift_register')
    
    for i = 1:sequence_length - 1
        seq(i) = prbs(i);
    end
    seq(sequence_length) = prbs(1);
    
else
    error('logical_prbs: pattern type not implemented.');
end

end