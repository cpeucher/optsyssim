function C = prod_mm(A,B)
% Product of time or frequency dependent 2x2 matrices
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements the multiplication between two 2x2 matrices that
% can depend on a parameter, typically time or frequency.
% This can be used, for instance, for Jones matrices calculations.
% The Jones matrices can represent time-varying polarization components,
% for instance, rotating elements or static polarization components.
% They can also be expressed as a function of frequency.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% C = prod_mm(A,B)
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% A                 matrix [complex 3D array]
%
%                       M is a 2x2xN array, where N is the number of
%                       samples: 
%                           M(:,:,kk) is the matrix for sample kk
%                           M(ii,jj,:) represents the evolution of the
%                           element ii, jj of the matrix.
%
% B                 matrix [complex 3D array]
%
%                       Same conventions are for matrix A
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% C                 matrix [complex 3D array]
%
%                       Same conventions are for input matrices
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

[na1,na2,na3]= size(A);
[nb1,nb2,nb3]= size(B);

if na1 ~= 2 || na2 ~= 2
    error('prod_mm: 2x2 matrix expected')
end

if nb1 ~= 2 || nb2 ~= 2
    error('prod_mm: 2x2 matrix expected')
end

if na3 ~= nb3
    error('prod_mm: number of samples of the 2 matrices is not consistent')
end

C = zeros(2,2,na3);
C(1,1,:) = squeeze(A(1,1,:)).'.*squeeze(B(1,1,:)).' + squeeze(A(1,2,:)).'.*squeeze(B(2,1,:)).';
C(1,2,:) = squeeze(A(1,1,:)).'.*squeeze(B(1,2,:)).' + squeeze(A(1,2,:)).'.*squeeze(B(2,2,:)).';
C(2,1,:) = squeeze(A(2,1,:)).'.*squeeze(B(1,1,:)).' + squeeze(A(2,2,:)).'.*squeeze(B(2,1,:)).';
C(2,2,:) = squeeze(A(2,1,:)).'.*squeeze(B(1,2,:)).' + squeeze(A(2,2,:)).'.*squeeze(B(2,2,:)).';


end