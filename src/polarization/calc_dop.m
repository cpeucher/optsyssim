function dop = calc_dop(S)
% Calculate degree of polarization from Stokes parameters
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the degree of polarization from the Stokes
% parameters.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% S = char_opt_stokes(sig);
% Sav = calc_stokes_average(S,nsamples);
% dop = calc_dop(Sav);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% S                 Stokes parameters [real matrix]
% 
%                       The input can either be a single Stokes vector,
%                       [S0; S1; S2; S3], or a set of Stokes vectors
%                       gathered in an 4xN matrix, where N is the number of
%                       samples:
%                           S(1,:) contains samples of Stokes parameter S0
%                           S(2,:) contains samples of Stokes parameter S1
%                           S(3,:) contains samples of Stokes parameter S2
%                           S(4,:) contains samples of Stokes parameter S3
% 
%                       The format of S is compatible with the output of
%                       the char_opt_stokes function, possibly averaged
%                       over time using the calc_stokes_average function.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% dop               degree of polarization [real vector]
% 
%                       dop is a line vector of dimension 1xN 
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

[m,~] = size(S);

if m ~= 4
    error('calc_dop: 4-element Stokes vector expected as input')
end

dop = sqrt(S(2,:).^2 + S(3,:).^2 + S(4,:).^2)./S(1,:);
% Degree of polarization

end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------