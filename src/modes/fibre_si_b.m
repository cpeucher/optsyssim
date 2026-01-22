function [b,nmodes] = fibre_si_b(mode,l,V,varargin)
% Calculate normalised propagation constant for circular-core step-index fibre
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function solves the dispersion equation of a given mode of a step
% index optical fibre. It returns the normalised propagation constant, 
% propagation constant and effective index for a given range of normalised
% frequencies V.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_fibre.a = 4.7e-6;
% params_fibre.n1 = 1.4628;
% params_fibre.n2 = 1.4600;
% lambda = 1550e-9;
% mode_type = 'LP';
% mode_l = 0;
% V = 2*pi*params_fibre.a*sqrt(params_fibre.n1^2 - params_fibre.n2^2)/lambda;
% [b, nmodes] = fibre_si_b(mode_type,mode_l,V,params_fibre);
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
% V                 normalised frequencies at which the dispersion equation
%                       will be solved [real vector]
%
% varargin          optional parameters: 
%                       optical fibre parameters [structure]
%
%                       varargin{1}.n2: refractive index of the cladding 
%                           [real scalar]
%
%                       varargin{1}.n1: refractive index of the core
%                           [real scalar]

% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% b                 normalised propagation constant
%
% nmodes            number of modes found for a given mode order (n for EH
%                   or HE, m for LP) at a given normalised frequency V.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


b_findroot = linspace(eps,1 - eps,1000);
% Normalised propagation constant axis for determination of the roots of
% the dispersion equation 


for ifreq = 1:length(V)    
    % Loop over normalised frequencies
    
    u_findroot = V(ifreq)*sqrt(1 - b_findroot);
    % Corresponding transverse wave number u
    
    f_root = fibre_si_dispersion_equation(mode,l,u_findroot,V(ifreq),varargin{:});
    % Calculate values of the eigenvalue equation at the b_findroot    
       
    cross = find(abs(abs(diff(sign(f_root)))-2)<1e-5);
    % Find indices where f_root changes sign

    
    nmodes(ifreq) = length(cross);
    % Number of roots of the function for a given normalised frequency
    
    % Now that we have the number of roots and their locations, we are
    % ready to solve the dispersion equation  
    
    if nmodes(ifreq) == 0
        b(1,ifreq) = NaN;
        % In case there are no roots at this particular we return NaN. 
        % But be careful not to count the size of the output array to 
        % determine the number of modes. Use nmodes instead.
    else  
        
        for iroot = 1:length(cross)
            % For each bracketed zero
            uroot(iroot,ifreq) = fzero(@(u)fibre_si_dispersion_equation(mode,l,u,V(ifreq),varargin{:}),[u_findroot(cross(length(cross)-iroot+1)) u_findroot(cross(length(cross)-iroot+1)+1)]);
            % Find root u of the dispersion equation
            b(iroot,ifreq) = (V(ifreq)^2 - uroot(iroot,ifreq)^2)/V(ifreq)^2;
            % Corresponding normalised propagation constant.
        end
        
    end
end    

maxroots = length(b);

for ifreq = 1:length(V)
    b(nmodes(ifreq)+1:maxroots,ifreq) = NaN;
end
% Clean up the arrays and save NaN below the cut-off of the modes. 
% Note the dirty way in which we proceeded, by "discovering" the number of
% roots for each frequency and adapting the size of the matrices each time
% a new mode is found...

end