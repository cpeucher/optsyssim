function sig = elec_sinusoidal(params)
% Electrical sinusoidal signal generation
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates a sinusoidal electrical signal.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_rf.frequency = symbol_rate;
% params_rf.phase = 0;
% params_rf.vpp = 1.0;
% params_rf.vdc = 0;
% sig = elec_sinusoidal(params_rf); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% params            electrical sinusoidal signal parameters [structure]
%
%                       params.frequency
%                           frequency of the signal, in Hz [real scalar]
%
%                       params.phase
%                           initial phase of the signal, in rad 
%                           [real scalar]
%
%                           For generating a cosine function: phase= 0;
%                           For generating a sine function:   phase= -pi/2;
%
%                       params.vpp
%                           peak-to-peak amplitude of the signal 
%                           [real scalar]
%
%                       params.vdc
%                           dc level of the signal [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               sinusoidal electrical signal [real vector]This is it.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% time_array        time samples, in s [real vector]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global time_array

sig = 0.5*params.vpp*cos(2*pi*params.frequency.*time_array + params.phase)+ params.vdc;

end
