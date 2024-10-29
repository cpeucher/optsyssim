function [freq_tf,magnitude,phase,phase_delay,freq_delay,group_delay] = elec_tf_analysis(freq,params,npoints,do_plot)
% Calculate and display electrical filter transfer function
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function analyzes the transfer function of an electrical filter by
% calculating and displaying its magnitude, phase, phase delay and group
% delay.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% freq              frequency range over which the transfer function will
%                       be analyzed, in Hz. 
%
%                       This can be a 2-point vector [fmin,fmax], or a full
%                       frequency vector [f0,f1,...,fN], in which case only
%                       f0 abd fN will be taken into account to define the
%                       interval of interest.
%
% params            structure containing parameters of the electrical 
%                       filter [structure]
%
% npoints           number of points over which the transfer function will
%                       be computed [integer]
%
% do_plot           switch to plot or not the transfer function
%
%                       do_plot = 0,1
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% freq_tf           frequency axis for magnitude, phase and phase delay,
%                       in Hz [real vector]
%
% magnitude         magnitude of the transfer funciton |h|, in linear scale
%                       [real vector]
%
%                       Take 20*log10() for value in dB.
%
% phase             phase of transfer function, in rad [real vector]
%
%                       It is here defined as arg(h)
%
% phase_delay       phase delay, in s [real vector]
%
% freq_delay        frequency axis for group delay, in Hz [real vector]
%
% group_delay       group delay, in s [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

freq_tf = linspace(freq(1),freq(end),npoints);
% Create frequency axis

tf = elec_tf_elpf(params,freq_tf);
% Calculate (complex) transfer function

magnitude = abs(tf);
phase = unwrap(angle(tf));

phase_delay = -phase./freq_tf/2/pi;
% Phase delay
% The minus sign stems from the fact that we define the phase as the
% argument of the transfer function (not minus the argument)

[freq_delay,group_delay] = num_diff(1,freq_tf,phase);
group_delay = -group_delay/2/pi;
% Group delay
% With our phase definition, the group delay is minus the derivative of the
% phase with respect to angular frequency.

if do_plot
    % Plot filter response.
    figure('Name','Transfer function analysis')
    subplot(2,2,1)
    plot(freq_tf,20*log10(magnitude),'b-')
    xlabel('frequency (Hz)')
    ylabel('|H| (dB)')
    title('magnitude')
    subplot(2,2,2)
    plot(freq_tf,phase,'b-')
    xlabel('frequency (Hz)')
    ylabel('\phi (rad)')
    title('phase')
    subplot(2,2,3)
    plot(freq_tf,phase_delay,'b-')
    xlabel('frequency (Hz)')
    ylabel('\tau_\phi (s)')
    title('phase delay')
    subplot(2,2,4)
    plot(freq_delay,group_delay,'b-')
    xlabel('frequency (Hz)')
    ylabel('\tau_g (s)')
    title('group delay')
end



end