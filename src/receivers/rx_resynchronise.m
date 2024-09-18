function [sig,imax] = rx_resynchronise(sig,seq,sig_type)
% Retiming of an electrical or optical signal based on intensity cross-correlation ("clock recovery")
% 
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function retimes a signal that has been delayed with respect to the
% initial binary sequence.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% sig_type = 'opt';%'elec';
% [sig,imax] = rx_resynchronise(sig,seq,sig_type);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical or electrical signal 
%                       [optical signal structure or real vector]
%
% seq               initial binary sequence [binary vector]
%
%                       The signal will be retimed to be in sync with this
%                       sequence.
%
% sig_type          type of the signal to retime [elec/opt]
%
%                       sig_type = 'opt';
%                       sig_type = 'elec'; 
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               retimed optical or electrical signal 
%                       [optical signal structure or real vector]
% 
% imax               number of samples the signal has been shifted in order
%                       to be resynchronised [integer scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

if strcmp(sig_type,'opt')
    nsamples = length(sig.x);
    
elseif strcmp (sig_type,'elec')
    nsamples = length(sig);
    
end
% Number of samples in the input signal

nsymbols = length(seq);
% Number of bits in the input pattern

nsps = nsamples/nsymbols;
% Number of samples per bit

refsig = seq(ones(1,nsps),:);
refsig = refsig(:).';
% Create a reference NRZ signal based on the input sequence data

if strcmp(sig_type,'opt')
    compare = abs(sig.x).^2 + abs(sig.y).^2;
    
elseif strcmp(sig_type,'elec')
    compare = sig;
    
end
% Determines which quantity is compared to the initial data, depending on
% whether the signal is optical or electrical.

[~,imax] = max(xcorr(compare,refsig));
% Calculate the delay shift that maximises the cross-correlation between 
% the signal and the reference signal.

if strcmp(sig_type,'opt')
    sig.x = circshift(sig.x,[0 -imax]);
    sig.y = circshift(sig.y,[0 -imax]);
    
elseif strcmp (sig_type,'elec')
    sig  =circshift(sig,[0 -imax]);
    
end
% Retime the signal

end