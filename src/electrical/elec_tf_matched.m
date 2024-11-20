function tf = elec_tf_matched(pulse,dt)
% Transfer function of matched electrical filter
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the transfer function of an electrical matched
% filter from the pulse shape.
% The transfer function is normalised so that the transmission at DC is
% equal to 1.
% In order to calculate the impulse response, the pulse is first retimed to
% the origin (t = 0), so that no (or little) delay is induced when applying 
% the matched filter.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% tf = elec_tf_matched(sig_pulse,dt);
% sig = elec_filter(sig,tf);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% pulse             signal pulse shape [complex vector]
%
% dt                time separation between samples, in s [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% tf                matched filter transfer function, in increasing 
%                       frequency order [complex vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% REMARKS:
% -------------------------------------------------------------------------
% The retiming operation should be pretty robust for standard "well-peaked"
% pulse shapes such as Gaussian, sech, sinc, raised-cosine, root
% raised-cosine, as well as flat-top pulse shapes.
% It has been thoroughly tested:
% 1. for peaked pulse when the peak maximum matches exactly one time
% sample.
% 2. for flat-top (rectangular) pulses with an odd number of samples, i.e.
% the center of gravity of the pulse exactly matches one time sample.
% It should work for other cases, possibly with a minor 1-sample delay or
% so due to lack of symmetry in the pulse definition.
% We did not take all cases into consideration. The function is likely to
% be used with "well-peaked" pulses such as root raised-cosine. 
% Furthermore, it is likely that some timing recovery operation will need
% to be performed at the receiver, which makes some possible 1-sample delay
% in the calculation of the matched filter impulse response irrelevant (it
% just looks nicer if we deal with this delay for back-to-back
% characterisation...)
% Therefore we did not spend more time investigating unlikely use cases.
% See the following test script for considerations on retiming the pulse 
% to the origin:
% \tests\manual\electrical\test_pulse_retiming.m
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

nsamples = length(pulse);
% Number of samples in the pulse

hr = fliplr(pulse(:).');
% Impulse response of the matched filter
% We first ensure the pulse shape is a line vector.
% Then we mirror it, i.e. the impulse response of the matched filter is
% hr(t) = K*p(-t), where p(t) is the pulse shape, and K is an arbitrary
% scaling parameter (see below).

K = 1/num_int1d_simpson(hr,dt);
% Normalisation factor so that the transfer function of the matched filter 
% at DC is equal to 1.

% We need to retime the pulse to the origin, so that no delay is induced
% when applying the matched filter.
% We first determine the "center of gravity" of the pulse by calculating
% the pulse 1st moment.
% This is more robust than simply looking for the max of the pulse, in case
% we have some flat-top pulse.

pulse_moment_sample = round(sum([1:nsamples].*hr)/sum(hr));
% Pulse first moment (digital).

hr = circshift(hr,-pulse_moment_sample + 1,2);
% Retime the impulse response so that it is centered at zero.

hr = K*hr;
% Normalised impulse response.

tf = num_ft(hr,dt);
% Transfer function of matched filter, in increasing frequency order.

end