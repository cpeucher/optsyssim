function [constellation,norm_es,norm_emax] = define_constellation(type,m)
% Define constellations of digital modulation formats
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function outputs a complex vector of length M, where M is the number
% of points in the constellation, containing the normalised complex values
% of all the points in the constellation.
%
% The constellation points are by defaut defined on a grid of spacing equal
% to 2, i.e. their real and imaginary parts take values in
% [... -9 -7 -5 -3 -1 1 3 5 7 9 ...].
% 
% The points in the constellation vector are indexed by the decimal value
% of the symbol + 1. For instance, for a 4-point constellation (2 bits per
% symbol):
% constellation(1) is the symbol associated with the transmission of
% decimal 0 / binary 00
% constellation(2) is the symbol associated with the transmission of
% decimal 1 / binary 01
% constellation(3) is the symbol associated with the transmission of
% decimal 2 / binary 10
% constellation(4) is the symbol associated with the transmission of
% decimal 3 / binary 11
%
% Normalisation factors are also provided so that the constellation has 
% an energy per symbol equal to 1, or a maximum energy equal to 1, i.e.
% that it fits within the circle of radius 1 
% The normalised constellation is obtained by 
% constellation/sqrt(norm_es) 
% constellation/sqrt(norm_emax)
%
% In order to obtain a normalised constellation, divide the output vector
% constellation by the normalisation factor of interest.
% Observe that this normalisation is not carried out within the function.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% [constellation,norm_es,norm_emax] = define_constellation(type,m);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% type              constellation to be generated [string]                 
%                       Currently defined type are:
%                       'bpsk'
%                       'pam4_natural'
%                       'pam4_gray'
%                       'pam8_natural'
%                       'pam16_natural'
%                       'qpsk_natural'
%                       'qpsk_gray'
%                       'psk8_natural'
%                       'psk8_gray'
%                       'qam8_rect_gray'
%                       'qam8_square_gray'
%                       'qam8_circ'
%                       'qam8_star'
%                       'psk16_natural'
%                       'psk16_gray'
%                       'qam16_natural'
%                       'qam16_quadrant'
%                       'qam16_gray'
%                       'qam16_differential'
%                       'qam32_rect_gray'
%                       'qam32_cross'
%                       'qam64_gray'
%
% m                 number of points in the constellation [integer]
%                       The number of bits per symbol is log2(m).
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% constellation     complex numbers representing the symbols in the
%                       constellation [complex vector]
%
% norm_es           energy per symbol of the constellation [real scalar]
%
% norm_emax         maximum energy of the constellation [real scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

symb = [0:m-1];
% Symbols indices for the given constellation
symb = dec2bin(symb);
% The symbols are now expressed in binary form

if rem(log2(m),2) == 0
    % Check whether the constellation is square
        msb = symb(:,1:log2(m)/2);
        % The first log2(M)/2 bits are defined as msb
        lsb = symb(:,log2(m)/2+1:log2(m));
        % The last log2(M)/2 bits are defined as lsb
        % OBS: valid for square constellations only!!
end


% ----
% Definition of constellations
% ----
switch type     
    
    case 'bpsk'
        
        constellation = [-1 1];
                
    case 'pam4_natural'
        
        constellation = [-3 -1 1 3];       
        
    case 'pam4_gray'
        
        constellation = [-3 -1 3 1];
        
    case 'pam8_natural'
        
        constellation = [-7 -5 -3 -1 1 3 5 7];
        
    case 'pam16_natural'
        
        constellation = [-15 -13 -11 -9 -7 -5 -3 -1 1 3 5 7 9 11 13 15];
        
    case 'pam32_natural'
        
        constellation = 2*(0:1:31) - 31;
        
    case 'pam64_natural'
        
        constellation = [-63:2:63];
        
    case 'qpsk_natural'
        
        mapx = [-1 1];
        mapy = [-1 1];
        
        constellation = mapx(not(bin2dec(msb)).*bin2dec(lsb) + bin2dec(msb).*not(bin2dec(lsb))+1) + 1i*mapy(bin2dec(msb)+1);
        % Create constellation     
        
    case 'qpsk_gray'
        
        mapx = [-1 1];
        mapy = [-1 1];
        
        constellation = mapx(bin2dec(lsb)+1) + 1i*mapy(bin2dec(msb)+1);
        
    case 'psk8_natural'
        
        constellation = [1 (1+1i)*sqrt(2)/2 1i (-1+1i)*sqrt(2)/2 -1 -(1+1i)*sqrt(2)/2 -1i (1-1i)*sqrt(2)/2];
        
    case 'psk8_gray'
        
        constellation = [1 (1+1i)*sqrt(2)/2 (1-1i)*sqrt(2)/2 1i -(1+1i)*sqrt(2)/2 -1 -1i -(1-1i)*sqrt(2)/2];
        
    case 'qam8_rect_gray'
        
        constellation = [-3+1i -3-1i -1+1i -1-1i 3+1i 3-1i 1+1i 1-1i];
        
    case 'qam8_square_gray'
        
        constellation = [-1-1i -1 -1i -1+1i 1 1+1i 1-1i 1i];          
        
    case 'qam8_circ'
        
        constellation = [1 exp(1i*6*pi/7) exp(1i*2*pi/7) exp(1i*4*pi/7) exp(1i*10*pi/7) exp(1i*8*pi/7) exp(1i*12*pi/7) 0];    
    
    
    case 'qam8_star'
        
        constellation = [-(1+sqrt(3)) -1+1i -1-1i 1i*(1+sqrt(3)) -1i*(1+sqrt(3)) 1+1i 1-1i 1+sqrt(3)];    
        
    case 'psk16_natural'
        
        constellation = [];        
        
    case 'psk16_gray'
        
        constellation = [];        
    
    case 'qam16_natural'
        % Natural mapping of QAM16 constellation

        mapx = [-3 -1 1 3];
        mapy = [-3 -1 1 3];
        % Mapping rule

        constellation = mapx(bin2dec(msb)+1) + 1i*mapy(bin2dec(lsb)+1);
        % Create constellation
        
        
    case 'qam16_quadrant'
        % Quadrant mapping (MSB indicates quadrant)
        
        mapmsb = [2+2i -2+2i -2-2i 2-2i];
        maplsb = [1+1i -1+1i -1-1i 1-1i];
        % Mapping rule

        constellation = mapmsb(bin2dec(msb)+1) + maplsb(bin2dec(lsb)+1);
        
        
    case 'qam16_gray'
        % Gray mapping of QAM16 modulation
        
        mapx=[-3 -1 3 1];
        mapy=[-3 -1 3 1];
        % Mapping rule

        constellation = mapx(bin2dec(msb)+1)+1i*mapy(bin2dec(lsb)+1);
        % Create constellation
        
%     case 'qam16_differential'
%         % Differential mapping. In progress. For the fun.       
%         
%         mapmsb = [2+2i -2+2i -2-2i 2-2i];
%         MAPLSB_Q1=[1+1i -1+1i -1-1i 1-1i];
%         
%         constellation = [];
%         
%         %iquad=bin2dec(msb);
%         
%         for iquad=0:3
%             % Loop over quadrant
%             
%             a = reshape(MAPLSB_Q1,[2 2]);
%             
%             a = rot90(a,iquad);
%             
%             a = reshape(a,[1 4]);
%             
%             constellation = [constellation mapmsb(iquad+1) + a(bin2dec(lsb) + 1)];
%         
%         end

    case 'qam32_rect_gray'
        
        constellation = [-7-3i -7-1i -7+3i -7+1i -5-3i -5-1i -5+3i -5+1i -1-3i -1-1i -1+3i -1+1i -3-3i -3-1i -3+3i -3+1i 7-3i 7-1i 7+3i 7+1i 5-3i 5-1i 5+3i 5+1i 1-3i 1-1i 1+3i 1+1i 3-3i 3-1i 3+3i 3+1i];

        
    case 'qam32_cross'
        
        constellation = [-3-5i -1-5i -3+5i -1+5i -5-3i -5-1i -5+3i -5+1i -1-3i -1-1i -1+3i -1+1i -3-3i -3-1i -3+3i -3+1i 3-5i 1-5i 3+5i 1+5i 5-3i 5-1i 5+3i 5+1i 1-3i 1-1i 1+3i 1+1i 3-3i 3-1i 3+3i 3+1i];
            
        
    case 'qam64_gray'
        
        % We'll clean up with msb lsb generation later...
        
        constellation = [-7-7i -7-5i -7-1i -7-3i -7+7i -7+5i -7+1i -7+3i -5-7i -5-5i -5-1i -5-3i -5+7i -5+5i -5+1i -5+3i -1-7i -1-5i -1-1i -1-3i -1+7i -1+5i -1+1i -1+3i -3-7i -3-5i -3-1i -3-3i -3+7i -3+5i -3+1i -3+3i 7-7i 7-5i 7-1i 7-3i 7+7i 7+5i 7+1i 7+3i 5-7i 5-5i 5-1i 5-3i 5+7i 5+5i 5+1i 5+3i 1-7i 1-5i 1-1i 1-3i 1+7i 1+5i 1+1i 1+3i 3-7i 3-5i 3-1i 3-3i 3+7i 3+5i 3+1i 3+3i];

        
    otherwise
        
        error('define_constellation: constellation not defined.');
        
end


% ----
% Calculation of normalisation factors
% ----
norm_es = sum(abs(constellation).^2)/m;
% Energy per symbol of the constellation
% Normalisation factor so that constellation/sqrt(normEs) is a unit energy
% constellation

% norm_es = sqrt(2*(m-1)/3);
% Theoretical normalisation factor for square constellation

norm_emax = max(abs(constellation).^2);
% Maximum energy of the constellation
% Normalisation factor so that the points of constellation/sqrt(normEmax) 
% have a maximum distance from the origin equal to 1

% norm_emax = (sqrt(2)*(sqrt(m)-1))^2.
% Theoretical normalisation factor for square constellation
% Assumes d = 2, which is the case of the constellations defined in this
% function

end