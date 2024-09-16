function sig = mod_iq(sig,sig_i1,sig_i2,sig_q1,sig_q2,v_ps,params)
% IQ modulator
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements a model of an optical IQ modulator made from two
% parallel Mach-Zehnder modulators and a phase shifter. 
% The model does not assume any particular operating point for the
% modulators and is therefore highly adaptable.
% A modulated signal can also be applied to the phase-shifter, making it
% possible to implement FSK modulation formats.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_iq.vpi_mzm = 1;
% params_iq.vpi_pm = 1;
% params_iq.bias_i1 = params_iq.vpi_mzm;
% params_iq.bias_i2 = 0;
% params_iq.bias_q1 = params_iq.vpi_mzm;
% params_iq.bias_q2 = 0;
% params_iq.split_i = 0.5;
% params_iq.split_q = 0.5;
% params_iq.loss_i = 0;
% params_iq.loss_q = 0;
% v_ps = params_iq.vpi_pm/2;
% sig_i1 = params_iq.vpi_mzm*(nrz_data_sig_i - 0.5);
% sig_i2 = -params_iq.vpi_mzmM*(nrz_data_sig_i - 0.5);
% sig_q1 = params_iq.vpi_mzm*(nrz_data_sig_q - 0.5);
% sig_q2 = -params_iq.vpi_mzm*(nrz_data_sig_q - 0.5);
% sig = mod_iq(sig,sig_i1,sig_i2,sig_q1,sig_q2,v_ps,params_iq);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
%                       
% sig_i1            voltage signal applied to the first arm of the in-phase
%                       Mach-Zehnder modulator [real vector]
%
% sig_i2            voltage signal applied to the second arm of the
%                       in-phase Mach-Zehnder modulator [real vector]
%
% sig_q1            voltage signal applied to the first arm of the 
%                       quadrature Mach-Zehnder modulator [real vector]
%
% sig_q2            voltage signal applied to the second arm of the 
%                       quadrature Mach-Zehnder modulator [real vector]
%
% v_ps              voltage applied to the phase shifter 
%                       [real scalar or real vector]
%                       If v_ps is a scalar, then a static phase shift is
%                       considered. Otherwise v_ps should be an electrical
%                       signal vector.
%
% params            IQ modulator parameters [structure]
%
%                       params.vpi_mzm
%                           half-wave voltage of the 2 Mach-Zehnder
%                           modulators [real scalar]
%
%                       params.vpi_pm 
%                           half-wave voltage of the phase shifter
%                           [real scalar]
%
%                       params.bias_i1
%                           bias applied to the first arm of the in-phase 
%                           Mach-Zehnder modulator [real scalar]
%
%                       params.bias_i2
%                           bias applied to the second arm of the in-phase 
%                           Mach-Zehnder modulator [real scalar]
%
%                       params.bias_q1
%                           bias applied to the first arm of the quadrature 
%                           Mach-Zehnder modulator [real scalar]
%
%                       params.bias_q2
%                           bias applied to the second arm of the 
%                           quadrature Mach-Zehnder modulator [real scalar]
%
%                       params.split_i
%                           power splitting ratio of the input and output
%                           Y-junctions of the in-phase Mach-Zehnder
%                           modulator [real scalar]
%
%                       params.split_q
%                           power splitting ratio of the input and output
%                           Y-junctions of the quadrature Mach-Zehnder
%                           modulator [real scalar]
%
%                       params.loss_i
%                           excess loss of the in-phase Mach-Zehnder 
%                           modulator, in dB [real scalar]
%
%                       params.loss_q
%                           excess loss of the quadrature Mach-Zehnder 
%                           modulator, in dB [real scalar]
%                           Use different values of
%                           params.loss_i and params.loss_q to
%                           simulate some imperfections in the input /
%                           output Y-junction of the global
%                           interferometric structure.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output optical signal [optical signal structure]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

[sig1,sig2] = opt_splitter_y_junction(sig);
% Input power splitter, assumbed to be ideal
sig1 = mod_mzm(sig1,sig_i1,sig_i2,params.bias_i1,params.bias_i2,params.vpi_mzm,params.split_i,params.split_i,params.loss_i);
% Data modulator for I quadrature.
sig2 = mod_mzm(sig2,sig_q1,sig_q2,params.bias_q1,params.bias_q2,params.vpi_mzm,params.split_q,params.split_q,params.loss_q);
% Data modulator for Q quadrature.
sig2 = opt_phase_shift(sig2,pi*v_ps/params.vpi_pm);
% pi/2 phase shift applied to the Q quadrature.
sig = opt_combiner_y_junction(sig1,sig2);
% Output power combiner.

end