function sig = mod_linear(sig,drive,extinction_ratio)
% Linear intensity modulator
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements an ideal linear optical modulator. It is assumed
% that both the -x and -y polarisations of the input signal are modulated 
% equally. 
% "Linear" refers to the power transfer of the modulator.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% extinction_ratio = 30;
% sig = mod_linear(sig,nrz_data_sig,extinction_ratio); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
%
% drive             driving electrical signal [real vector]
%                       Ideally the driving signal should take values 
%                       between 0 and 1, in which case the operation of the 
%                       modulator is quite straightforward.
%                       Otherwise things may get more complicated since
%                       one may have an amplifying modulator, or possibly 
%                       truncation of the modulating signal to ensure that 
%                       the argument of the sqrt in the field transfer 
%                       function remains positive. 
%                       (Some) warning messages are issued.
%
% extinction_ratio  desired extinction ratio of the modulated signal,
%                       in dB [real scalar]
%                       Currently the value "+Inf" is not supported for the 
%                       extinction ratio. But a value of 100 or even 
%                       1000 dB should be sufficient for most purposes...
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
% REMARKS:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% TO DO:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% CREDITS:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% AUTHOR:
% -------------------------------------------------------------------------
% Christophe Peucheret (christophe.peucheret@univ-rennes1.fr)
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

pin = abs(sig.x).^2+abs(sig.y).^2;
% Power of the input signal.
pin_peak = max(pin);
% Peak power of the input signal.

extinction_ratio = 10^(extinction_ratio/10);
% Convert extinction ratio from dB to linear scale.

d0 = min(drive);
if d0 ~= 0
    Message{1} = ['mod_linear (warning): Minimum of driving signal= ' num2str(d0,'%3.2e') ' (~= 0)' ' Delta= ' num2str(abs(d0),'%3.2e') '.'];
else
    Message{1} = ['mode_linear: Minimum of driving signal= ' num2str(d0)];
end

d1 = max(drive);
if d0 ~= 0
    Message{2} = ['mod_linear (warning): Maximum of driving signal= ' num2str(d1,'%3.2e') ' (~= 1)' ' Delta= ' num2str(abs(d1-1),'%3.2e') '.'];
else
    Message{2} = ['mod_linear: Maximum of driving signal= ' num2str(d1)];
end

modindex = (extinction_ratio-1) / (d1 - 1 + extinction_ratio - extinction_ratio*d0);
% Calculate modulation index from the required extinction ratio and minimum
% and maximum values of the driving signal.

modulation=(1 - modindex) + modindex*drive;
% Calculate modulation signal for the required modulation index.

test = (modulation<0);

if sum(test)~=0
    Message{3} = 'mod_linear (warning): Modulation signal is negative. Will be truncated to zero.';
    modulation(test) = 0;
else
    Message{3} = 'mod_linear: Modulation signal is always positive (OK).';
end
% Truncate the modulation signal in case it takes negative values.


sig.x = sig.x.*sqrt(modulation);
sig.y = sig.y.*sqrt(modulation);
% Apply modulating signal to the -x and -y polarisations of the input
% field.

pout = abs(sig.x).^2+abs(sig.y).^2;
% Power of the output signal.
pout_peak = max(pout);
% Peak power of the output signal.

if pout_peak > pin_peak
    Message{4} = ['mod_linear (warning): Output peak power (' num2str(pout_peak/1.0e-3) ' mW) larger than input peak power (' num2str(pin_peak/1.0e-3) ' mW).'];
else
    Message{4} = 'mod_linear: Output peak power at most equal to input peak power (OK).';
end

pout_min = min(pout);

er_calc = pout_peak/pout_min;
% Calculate the extinction ratio at the output of the modulator. Observe
% that this definition is suitable only with "well behaved" driving signal.

end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------