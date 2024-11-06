function varargout = conv_loss_lin_log(loss,varargin)
% Conversion of loss per unit length from dB/km to m^-1 and calculation of effective length
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function converts loss expressed in the usual engineering unit of
% dB/km to the SI unit of m^-1 (sometimes refered to as nepers).
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% loss = 0.8;       % Loss in dB/km
% length = 50;      % Waveguide length, in m.
% [alpha,leff] = conv_loss_lin_log(loss,length);
% fprintf('\n\n%s\t\t\t\t%3.2e\t\t%s\n','Loss',loss,'dB/km');
% fprintf('%s\t\t\t\t%3.2e\t\t%s\n','Loss',alpha,'1/m');
% fprintf('%s\t\t\t\t%3.2e\t\t%s\n','Length',length,'dB/km');
% fprintf('%s\t%3.2e\t\t%s\n','Effective length',leff,'m');
% fprintf('%s\t\t\t\t%3.2e\t\t%s\n','1/alpha',1/alpha,'m');
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% loss              loss, in dB/km [real vector]
%                       The number should be positive for a lossy medium,
%                        e.g. 0.2 dB/km for standard fibre.
%
% varargin          fibre length, in m [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% alpha             loss per unit length, in m^-1 [real vector]
%
% leff              fibre effective length, in m [real vector] 
%                       Only when the fibre length is provided as an
%                       optional input argument.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

alpha = loss*log(10)/10000;
% Convert loss to 1/m

optargin = size(varargin,2);
% Number of optional arguments
if optargin == 0
    varargout{1} = alpha;
elseif optargin == 1
    if alpha == 0
        leff = varargin{1};
    else
        leff = (1 - exp(-alpha*varargin{1}))./alpha;
        % Effective length in m
    end
    varargout{1} = alpha;
    varargout{2} = leff;
else
    error('conv_loss_lin_log: too many input arguments.');
end

end