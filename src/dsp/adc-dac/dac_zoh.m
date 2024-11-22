function xa = dac_zoh(x,ndac)
% Basic zero-order hold interpolation
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements a basic zero-order hold interpolation by
% repeating each input sample by the the specified oversampling factor.
% Other options are implemented elsewhere, for instance up-sample using 
% dsp_upsample.m followed by electrical low pass-filtering
% with a sinc filter (type = 'zero_order_hold' in elec_elpf.m).
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% sig = dac_zoh(samples,ndac);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% x                 input samples [complex vector]
%
% ndac              oversampling factor [integer]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% xa                sampled analog waveform [complex vector]
%                       The length of the vector is length(x)*ndac.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

xa = x(ones(1,ndac),:);
xa = xa(:).';


end