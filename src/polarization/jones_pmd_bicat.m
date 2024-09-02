function T = jones_pmd_bicat(nsections,pmd,h,alpha,phi,omega)
% Jones matrix of concatenation of birefringent elements and rotators
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements ...
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% nsection          number of birefringent sections [integer scalar]
%                       This is it.
%
% alpha             polarisation rotation angle of each section 
%                       [real vector]
% 
% phi               phase shift [real vector]
%

% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output optical signal [optical signal structure]
%                       This is it.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
%
%
% CONSTANT              essential physical constants [structure]
%
% time_array            time samples, in s [real vector]
%
% dt                    time samples separation, in s [real scalar]
%
% frequency_array       relative frequency samples, in Hz [real vector]
%
% reference_frequency   reference frequency, in Hz [real scalar]
%
% df                    frequency samples separation, in Hz [real scalar]
%
% space_grid            2D space grid in the transverse plane [structure]
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
% Christophe Peucheret (christophe.peucheret@univ-rennes.fr)
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
nsamples = length(omega);


T = zeros(2,2,nsamples);
T(1,1,:) = ones(1,nsamples);
T(2,2,:) = ones(1,nsamples);
% Initialisation of the Jones matrix of the link (identity matrix)


for iseg = 1:nsections
    % Loop over birefringent fibre segments
    
    B = zeros(2,2,nsamples);
    % Initialise birefringent element matrix of the section.
    R = zeros(2,2,nsamples);
    % Initialise rotation matrix of the section.
    
    R(1,1,:) = cos(alpha(iseg));
    R(1,2,:) = sin(alpha(iseg));
    R(2,1,:) = -sin(alpha(iseg));
    R(2,2,:) = cos(alpha(iseg));
    
    B(1,1,:) = exp(1i*(sqrt(3*pi/8)*pmd*omega*sqrt(h(iseg))/2 + phi(iseg)));
    B(2,2,:) = conj(B(1,1,:));
    
    Tspan = prod_mm(B,R);
    % Jones matrix of the segement 
    
    T = prod_mm(Tspan,T);  
    % Update Jones matrix of the link.
    
end
% End of loop over birefringent fiber segments





end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------