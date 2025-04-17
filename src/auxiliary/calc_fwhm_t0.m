function out = calc_fwhm_t0(pulse_type,pulse_order,mode,in)
% Conversion between different definitions of pulse widths
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function converts pulse widths from full-width at half-maximum to
% usual T0 for classic Gaussian and sech pulse shapes and various pulse 
% orders.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% pulse_type = 'sech'; %'gaussian';
% pulse_order = 1;
% pulse_width_fwhm = 10e-12;
% pulse_width_t0 = calc_fwhm_t0(pulse_type,pulse_order,'from_fwhm',pulse_width_fwhm);
% pulse_width_fwhm = calc_fwhm_t0(pulse_type,pulse_order,'to_fwhm',pulse_width_t0);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% pulse_type        pulse type [string]
%                       pulse_type = 'gaussian';
%                           u(t) = exp(-0.5*(t/t0)^2)
%                       pulse_type = 'sech';
%                           u(t) = sech(t/t0)
%
% pulse_order       pulse order [integer]
%
% mode              converts from or to the FWHM [string]
%                       mode = 'to_fwhm';
%                       mode = 'from_fwhm';
%
% in                value of the pulse width to convert [real scalar]
%                       T0 pulse width   if mode = 'to_fwhm'.
%                       FWHM pulse width if mode = 'from_fwhm'.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% out               converted pulse width, in the same unit as the input
%                       pulse width [real scalar] 
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

if strcmp(pulse_type,'gaussian')
    
    factor = 2*(log(2))^(1/(2*pulse_order));
    
    if strcmp(mode,'to_fwhm')
        out = factor*in;        
    elseif strcmp(mode,'from_fwhm')
        out = in/factor;        
    end   
    
elseif strcmp(pulse_type,'sech')
    
    if pulse_order == 1        
        factor = 2*log(1+sqrt(2));
        % calculate t0 from the expected FWHM
        % see e.g. Agrawal, Nonlinear Fiber Optics, 3rd edition, Academic
        % Press, pp. 71-72, 2001.        
        if strcmp(mode,'to_fwhm')
            out = factor*in;
        elseif strcmp(mode,'from_fwhm')
            out = in/factor;
        end        
    else
        error('calc_fwhm_t0: order currently limited to 1 for sech pulses.');
    end
    
else
    disp('calc_fwhm_t0: pulse type not implemented.');
end

%--------------------------------------------------------------------------
% End of function
%--------------------------------------------------------------------------