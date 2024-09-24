function [sig21,sig22] = opt_mzdi(sig11,sig12,params)
% Mach-Zehnder delay interferometer
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements a delay interferometer to be used for e.g. the
% optical demodulation of DnPSK signals. A Mach-Zehnder structure is
% implemented. In the model, a delay is placed in the upper arm of the
% interferometer, while the interferometer can be tuned by applying a phase
% shifting element to its lower arm. The interferometer is assumed to be
% strictly polarisation insensitive.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_mzdi.delay_target = 1/symbol_rate;
% params_mzdi.coupling_ratio_in = 0.5;
% params_mzdi.coupling_ratio_out = 0.5;
% params_mzdi.mode = 'tuned';%'general';
% params_mzdi.phase_shift = 0;
% [sig21,sig22] = opt_mzdi(sig,opt_nosig,params_mzdi); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig11             optical signal applied to the upper input port
%                       [optical signal structure]
%
% sig12             optical signal applied to the lower input port
%                       [optical signal structure]
%
% params            MZDI parameters [structure]
%
%                       params.delay_target      
%                           target interferometer delay, in s [real scalar]
% 
%                       params.coupling_ratio_in
%                           input coupler coupling ratio, in linear units
%                           [real scalar]
% 
%                       params.coupling_ratio_out
%                           output coupler coupling ratio, in linear units
%                           [real scalar]
%
%                       params.mode              
%                           specifies whether the implemented 
%                           interferometer is  made tuned to the signal 
%                           frequency, or if a phase shift corresponding to
%                           the delay is also applied to the upper arm, in 
%                           which case tuning requires a proper phase 
%                           tuning of the lower arm [string]
%
%                           mode = 'tuned';
%                           mode = 'general';
%
%                       params.phase_shift       
%                           value of the phase shift applied to the lower 
%                           arm [real scalar]
%
%                           Observe that this phase shift is applied even 
%                           when the delay interferometer is operated in
%                           the 'tuned' mode, thus allowing the 
%                           demodulation of formats such as DQPSK, D8PSK,
%                           etc in this mode, for instance.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig21             optical signal at the upper output port
%                       [optical signal structure]
%
% sig22             optical signal at the lower output port
%                       [optical signal structure]
%
% -------------------------------------------------------------------------
% REMARKS:
% -------------------------------------------------------------------------
% The params structure is updated with:
%
%                       params.delay_actual      
%                           actual delay that has been applied, 
%                           corresponding to an integer number of samples
%                           [real scalar]
%
%                           delay_actual is an integer multiple of dt
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

if strcmp(mode,'tuned')
    include_delay_phase_shift = 0;
elseif strcmp(mode,'general')
    include_delay_phase_shift = 1;
else
    disp('opt_mzdi: phase shift option not implemented.')
end

[sigup,siglow] = opt_coupler_2x2(sig11,sig12,'lin',params.coupling_ratio_in);
% Input 3 dB coupler

[sigup,params.delay_actual] = opt_delay(sigup,include_delay_phase_shift,params.delay_target);
% Optical delay in applied to the upper arm of the interferometer

siglow = opt_phase_shift(siglow,params.phase_shift);
% Phase shift applied to the lower arm of the interferometer

[sig21,sig22] = opt_coupler_2x2(sigup,siglow,'lin',params.coupling_ratio_out);
% Output 3 dB coupler.

end