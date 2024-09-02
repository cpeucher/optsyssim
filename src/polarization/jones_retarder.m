function R = jones_retarder(retardation)
% Jones matrix of a wave retarder
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements the Jones matrix of a wave retarder whose 
% slow axis coincides with the -x axis of the canonical basis (x,y) and
% whose retardation is specified by the parameter retardation.
% In particular:
% for retardation = 2*pi   ->   wave plate
% for retardation = pi     ->   half wave plate
% for retardation = pi/2   ->   quarter wave plate
% The retardation can be time-dependent.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% retardation       phase retardation, in radians [real vector]
%
%                       In a real implementation of a wave plate,this would 
%                       correspond to 
%                       retardation = 2*pi/lambda (n_s -n_f)*d 
%                       where n_s and n_f are the refractive indices
%                       associated with propagation along the slow and fast
%                       axis, respectively, and d is the thickness
%                       of the plate.
% 
%                       For a static retarder, retardation is a scalar.
%   
%                       For a dynamic retarder, retardation is a line 
%                       vector of dimension 1xN, where N is the number of
%                       samples.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% R                 Jones matrix of the wave retarder [2x2xN complex array]
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
% Christophe Peucheret (christophe.peucheret@univ-rennes.fr)
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

retardation = retardation(:).';
% Ensure retardation is a line vector.

% R = [exp(1i*retardation/2)         0;
%                  0              exp(-1i*retardation/2)];
% Jones matrix of wave retarder. 

R(1,1,:) = exp(1i*retardation/2);
R(1,2,:) = zeros(1,length(retardation));
R(2,1,:) = zeros(1,length(retardation));
R(2,2,:) = exp(-1i*retardation/2);

end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------