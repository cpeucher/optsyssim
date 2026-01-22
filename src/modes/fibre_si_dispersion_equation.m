function f = fibre_si_dispersion_equation(mode,l,u,V,varargin)
% Dispersion equations for the modes in a perfectly circular step-index
% optical fibre
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements dispersion equations for the modes in a step-
% index fibre with perfectly circular core geometry.
% Typically, this function will be called to calculate the
% propagation constants of the modes for a given normalised frequency.
% See e.g. K. Okamoto, Fundamentals of Optical Waveguides, 2nd ed.
% (Elsevier, 2006)
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% mode              mode family for which the dispersion equation should be 
%                       solved [string]
%
%                       mode = 'TE','TM','HE','EH','LP'
%           
%                       For the time being, exact dispersion equations are
%                       implemented for TE and TM modes (in the latter
%                       case, this requires the refractive indices of the 
%                       core and the cladding to be passed as an optional
%                       parameter) while approximations under the weakly-
%                       guiding approximation are implemented for HE and EH 
%                       modes.
%
% l                 azimutal mode index (integer scalar)
%
%                       In case mode = 'HE' or 'EH', the  argument l >= 1
%                       should be provided to calculate the HElm or EHlm
%                       modes.
%
%                       In case mode = 'TE' or 'TM', the  argument l = 0
%                       should be provided to calculate the TE0m or TM0m
%                       modes.
%
%                       In case mode = 'LP', the  argument l >= 0
%                       should be provided to calculate the LPlm modes 
%
%                       (l = azimutal mode order; m = radial mode order)
%
% u                 transverse wave number [real vector]
%                   
%                       Argument of the dispersion equation. 
%
% V                 normalised frequency [real scalar]
%
% varargin          optional parameters: 
%                       optical fibre parameters [structure]
%
%                       varargin{1}.n2: refractive index of the cladding 
%                           [real scalar]
%
%                       varargin{1}.n1: refractive index of the core
%                           [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% f                 value of the dispersion equation for the transverse 
%                       wave number u
% 
%                       The roots u of f(u) will determine the propagation 
%                       constants of the modes.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


if rem(l,1) ~= 0
        error('fibre_si_dispersion_equation: the azimutal mode number should be an integer.');
end

if (strcmp(mode,'HE') | strcmp(mode,'EH')) & (l < 1)
    error('fibre_si_dispersion_equation: the azimutal mode number l should be strictly positive l >= 1 for EH and HE modes.');
end

if (strcmp(mode,'TE') | strcmp(mode,'TM')) & (l ~= 0)
    error('fibre_si_dispersion_equation: the l parameter should be equal to 0 for TE and TM modes.');
end

if strcmp(mode,'LP') & (l < 0)
    error('fibre_si_dispersion_equation: the azimutal mode number l should be positive l >= 0 for LP modes.');
end


noptargin = size(varargin,2);
% Number of optional arguments

if noptargin > 1
    error('fibre_si_dispersion_equation: the number of optional arguments is limited to 1. It should be a structure containing fibre parameters for TM modes');
    % Only one optional argument is allowed
end

if strcmp(mode,'TM') 
    if noptargin == 0
    error('fibre_si_dispersion_equation: fibre parameters are required for TM modes.');
    else 
    params = varargin{1};
    end
end


if u < 0
    error('fibre_si_dispersion_equation: the transverse wave number u should be positive');
end
% Test input arguments


w = sqrt(V.^2 - u.^2);
% Calculate the transverse wave number w from the transverse wave number u
% and the normalised frequency V


switch mode
    % Switch over mode type

    case 'TE'
            f = besselj(1,u).*besselk(0,w).*w + besselj(0,u).*besselk(1,w).*u;
            % Okamoto, eq. (3.19)
            % The roots u of the function f should be found to calculate
            % the propagation constants of the TE0m modes.
            % Note that this is an exact solution. 
    
    case 'TM'
            f = (params.n2/params.n1)^2*besselk(1,w).*besselj(0,u).*u + besselk(0,w).*besselj(1,u).*w;
            % Okamoto, eq. (3.28)
            % The roots u of the function f should be found to calculate 
            % the propagation constants of the TM0m modes.
            % Note that this is an exact solution. 
            % Under the weakly guiding approximation, this solution reduces
            % to the TE equation below.         
        
    case 'HE'
            f = besselj(l - 1,u).*besselk(l -2,w).*w + besselj(l - 2,u).*besselk(l - 1,w).*u;
            % Okamoto, eq. in Table 3.1
            % The roots u of the function f should be found to calculate 
            % the propagation constants of the HElm (l >= 1) modes.
            % The case HE1l does not need to be considered separately as in
            % Okamoto since besselj(-1,x) = - besselj(1,x) and
            % besselk(-1,x) = besselk(1,x)
            % Note that this is an approximated solution under the weakly guiding
            % approximation.         
        
    case 'EH'
            f = besselj(l + 1,u).*besselk(l,w).*w + besselj(l,u).*besselk(l + 1,w).*u;
            % Okamoto, eq. in Table 3.1
            % The roots u of the function f should be found to calculate the
            % propagation constants of the EHlm (l >= 1) modes.
            % Note that this is an approximated solution under the weakly guiding
            % approximation.           
                    
    case 'LP'
            f = besselj(l,u).*besselk(l - 1,w).*w + besselj(l - 1,u).*besselk(l,w).*u;
            % Okamoto, eq. (3.74)
            % The roots u of the function f should be found to calculate 
            % the propagation constants of the LPlm (l >= 0) modes.
            % Note that this mode definition is only valid under the weakly 
            % guiding approximation.    
        
    otherwise
            error('fibre_si_dispersion_equation: mode type does not exist.');       
        
end
% End of switch over mode type

end