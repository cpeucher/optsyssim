function B = prod_mv(M,A)
% Product of time or frequency dependent 2-element vector and 2x2 matrix
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements the multiplication between a 2x2 matrix and a
% 2x1 vector. Both vector and matrix can depend on a parameter, typically
% time or frequency.
% This can be used for instance for Jones vector calculations: 
% - The Jones vector can consist of samples of a time-varying signal. 
% - The Jones matrix can represent a time-varying polarization component,
% for instance a rotating element) or a static polarization component.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% B = prod_mv(M,A)
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% M                 matrix [complex 3D array]
%
%                       M is a 2x2xN array, where N is the number of
%                       samples: 
%                           M(:,:,kk) is the matrix for sample kk
%                           M(ii,jj,:) represents the evolution of the
%                           element ii, jj of the matrix.
%
% A                 input vector [complex 2D array]
%
%                       A(:,kk) is the input vector for sample kk
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% B                 output vector [complex 2D array]
%
%                       Same conventions as for the input vector
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

[nlv,ncv] = size(A);
[nm1,nm2,nm3]= size(M);

if nlv ~= 2
    error('prod_mv: 2-dimensional vector expected')
end

if nm1 ~= 2 || nm2 ~= 2
    error('prod_mv: 2x2 matrix expected')
end

if nm3 > 1
    % The matrix is parameter-dependent.
    if ncv ~= nm3
        error('prod_mv: inconsistency between number of elements of vector and matrix')
    end
end
    
B = [squeeze(M(1,1,:)).'.*A(1,:) + squeeze(M(1,2,:)).'.*A(2,:);
     squeeze(M(2,1,:)).'.*A(1,:) + squeeze(M(2,2,:)).'.*A(2,:)];

end