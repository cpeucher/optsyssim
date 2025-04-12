function [pse,status] = sep_awgn_reference(esn0_db,mod,m)
% Reference SEP versus Es/N0 curves over AWGN for various constellations 
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function returns theoretical symbol-error probabilities as a 
% function of signal-to-noise ratio, specified as Es/N0, for some standard
% constellations over an additive white Gaussian noise (AWGN) channel.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% [pse, 'status'] = sep_awgn_reference(esn0_db,'qpsk',m); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% esn0_db           energy per symbol devided by noise spectral density,
%                       in dB [real vector]
%
% mod               modulation format [string]
%
%                       mod = 'bpsk'
%                       mod = 'qpsk'
%                       mod = 'pam_antipodal'
%                       mod = 'pam'
%                       mod = 'qam_square' 
%
% m                 number of symbols in the considered constellation
%                       [integer scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% pse               symbol-error probability [real vector]
%
% status            specifies the exact or approximate nature of the
%                       returned symbol error probability [string]
%
%                       status = 'exact'
%                       status = 'approximation'
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% TO DO:
% -------------------------------------------------------------------------
% Check my derivations for exact / approximation status.
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

esn0 = 10.^(esn0_db/10);
% Es/N0 in linear scale

switch mod    
    
    case 'bpsk'
        
        pse = func_q(sqrt(2*esn0));
        status = 'exact';
        
        
    case 'qpsk'
        
        qq = func_q(sqrt(esn0));
        
        pse = 2*qq - qq.*qq;
        status = 'exact';
        
        
    case 'pam_antipodal'
        
        pse = 2*(m - 1)/m * func_q(sqrt(6*esn0/(m^2 - 1)));
        status = '-';
        
    case 'pam'
        
        pse = 2*(m - 1)/m * func_q(sqrt(3*esn0/(m - 1)/(2*m-1)));
        status = '-';
        
        
    case 'qam_square'
        
        %         qq = func_q(sqrt(3*EsN0/(M - 1)));
        %
        %         pse = 4*(sqrt(M) - 1)*qq/sqrt(M) - 4*((sqrt(M)-1)/sqrt(M))^2*qq.*qq;
        
        pse = 1 - (1 - (sqrt(m) - 1)/sqrt(m)*erfc(sqrt(3*esn0/(m - 1)/2))).^2;
        status = '-';
        
    otherwise
        
        error('sep_awgn_reference: modulation format not defined.');
        
end 
% End of switch over modulation type
        
        
end