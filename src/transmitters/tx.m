function sig = tx(params)
% General optical transmitter
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function is to be used as a calling function to various standard
% optical transmitters.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_tx.type = 'ook_nrz';%'ook_rz33', 'ook_rz50', 'ook_rz67',
% %                       'ook_db_nrz_delay_add', 'ook_db_nrz_low_pass'
% %                       'dpsk_nrz_mzm', 'dpsk_nrz_pm'
% %                       'dpsk_rz33_mzm', 'dpsk_rz50_mzm', 'dpsk_rz67_mzm'
% %                       'qpsk_nrz_pm_pm', 'qpsk_nrz_mzm_pm', 'qpsk_nrz_iq'
% %                       'pam4_nrz_mzm'
% params_tx.emission_frequency = reference_frequency;
% params_tx.power = 1.0e-3;
% params_tx.linewidth = 0;
% params_tx.bit_pattern = bit_pattern;
% params_tx.bit_pattern_1 = bit_pattern_1;
% params_tx.bit_pattern_2 = bit_pattern_2;
% params_tx.rise_time = 1/symbol_rate/4;
% sig = tx(params_tx);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% params            transmitter parameters [structure]
%
%                       params.type
%                           transmitter type [string]
%                           params.type = 'ook_nrz'
%                           params.type = 'ook_rz33'
%                           params.type = 'ook_rz50'
%                           params.type = 'ook_rz67'
%                           params.type = 'ook_db_nrz_delay_add'
%                           params.type = 'ook_db_nrz_low_pass'
%                           params.type = 'dpsk_nrz_mzm'
%                           params.type = 'dpsk_nrz_pm'
%                           params.type = 'dpsk_rz33_mzm'
%                           params.type = 'dpsk_rz50_mzm'
%                           params.type = 'dpsk_rz67_mzm'
%                           params.type = 'qpsk_nrz_pm_pm'
%                           params.type = 'qpsk_nrz_mzm_pm'
%                           params.type = 'qpsk_nrz_iq'
%                           params.type = 'pam4_nrz_mzm'   
%
%                       params.emission_frequency
%                           emission frequency of the transmitter, in Hz
%                           [real scalar]
%
%                       params.bit_pattern
%                           bit pattern to be modulated, for binary 
%                           transmitters [binary vector]
%
%                       params.power
%                           power of the laser, in W [real scalar]
%
%                       params.linewidth
%                           laser linewidth, in Hz [real scalar]
%
%                       params.bit_pattern_1
%                           bit pattern to be modulated, for quaternary
%                           transmitters [binary vector]
%
%                       params.bit_pattern_2
%                           bit pattern to be modulated, for quaternary 
%                           transmitters [binary vector]
%
%                       params.bit_pattern
%                           structure containg the 4 bits patterns for 
%                           16-ary transmitter (16QAM) [structure]
%
%                       params.rise_time
%                           rise time of the electrical signal driving the 
%                           modulator, in s [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               modulated optical signal [optical signal structure]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

switch params.type
    % Switch over transmitter type (modulation format)
    
    case 'ook_nrz'
        
        sig = tx_ook_nrz(params);
        
    case 'ook_rz33'
        
         sig = tx_ook_rz33(params);
        
    case 'ook_rz50'
        
         sig = tx_ook_rz50(params);
         
    case 'ook_rz67'
        
         sig = tx_ook_rz67(params);
         
    case 'ook_db_nrz_delay_add'
        
        sig = tx_duobinary_nrz_delay_add(params);        
        
    case 'ook_db_nrz_low_pass'
        
        sig = tx_duobinary_nrz_low_pass(params);
        
    case 'dpsk_nrz_mzm'
        
        sig = tx_dpsk_nrz_mzm(params);
        
    case 'dpsk_nrz_pm'
        
        sig = tx_dpsk_nrz_pm(params);
        
    case 'dpsk_rz33_mzm'
        
        sig = tx_dpsk_rz33_mzm(params);
        
    case 'dpsk_rz50_mzm'
        
        sig = tx_dpsk_rz50_mzm(params);
        
    case 'dpsk_rz67_mzm'
        
        sig = tx_dpsk_rz67_mzm(params);
        
    case 'qpsk_nrz_pm_pm'
        
        sig = tx_qpsk_nrz_pm_pm(params);
        
    case 'qpsk_nrz_mzm_pm'
        
        sig = tx_qpsk_nrz_mzm_pm(params);
        
    case 'qpsk_nrz_iq'
        
        sig = tx_qpsk_nrz_iq(params);
        
    case 'pam4_nrz_mzm'
        
        sig = tx_pam4_mzm(params);
        
    otherwise
        
        error('tx: transmitter type not defined.');
        
end
% End of switch over modulation format.

end