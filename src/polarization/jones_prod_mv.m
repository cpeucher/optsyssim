function B = jones_prod_mv(M,A)
% Product of a Jones vector by a Jones matrix (possibly time-dependent)
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements the multiplication between a Jones matrix and a
% Jones vector. Both vector and matrix can be time dependent. 
% The Jones vector can consist of samples of a time-varying signal. 
% The Jones matrix can represent a time varying polarization component (for
% instance a rotating element) or a static polarization component.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% B = jones_prod_mv(M,A)
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% M                 Jones matrix [complex 3D array]
%
%                       For a static polarization component:
%                           M is a 2x2 Jones matrix [complex matrix]
%
%                       For a dynamic polarization component:
%                           M is a 2x2xN array, where N is the number of
%                           time samples.
%                           M(:,:,kk) is the Jones matrix at instant kk
%                           M(ii,jj,:) represents the time evolution of the
%                           element ii, jj of the Jones matrix.
%                       M is represented on the canonical basis x, y
%
% A                 input Jones vector [complex 2D array]
%
%                       A(:,kk) is the input Jones vector at instant kk
%                       A is represented on the canonical basis x, y
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% B                 output Jones vector [complex 2D array]
%
%                       Same conventions as for the input Jones vector
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

[nlv,ncv] = size(A);
[nm1,nm2,nm3]= size(M);

if nlv ~= 2
    error('jones_prod_mv: 2-dimensional vector expected')
end

if nm1 ~= 2 || nm2 ~= 2
    error('jones_prod_mv: 2x2 matrix expected')
end

if nm3 > 1
    % The Jones matrix is time-dependent.
    if ncv ~= nm3
        error('jones_prod_mv: inconsistency between number of elements of vector and matrix')
    end
end
    
B = [squeeze(M(1,1,:)).'.*A(1,:) + squeeze(M(1,2,:)).'.*A(2,:);
     squeeze(M(2,1,:)).'.*A(1,:) + squeeze(M(2,2,:)).'.*A(2,:)];

end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------