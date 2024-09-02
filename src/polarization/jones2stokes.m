function S = jones2stokes(E)
% Conversion from Jones to Stokes vector
%
% Parameters
% ==========
% E : complex 2xN matrix
%     Matrix of Jones vectors
%     Each column of E is a 2-element Jones vector
%
% Returns
% =======
% S : real 4xN matrix
%     Matrix of Stokes vectors (including S0)
%     Each column of S is a 4-element Stokes vector


[~,nvec] = size(E);

S = zeros(4,nvec);

for ii = 1:nvec
    
    S(1,ii) = E(1,ii)*conj(E(1,ii)) + E(2,ii)*conj(E(2,ii));
    S(2,ii) = E(1,ii)*conj(E(1,ii)) - E(2,ii)*conj(E(2,ii));
    S(3,ii) = E(1,ii)*conj(E(2,ii)) + E(2,ii)*conj(E(1,ii));
    S(4,ii) = 1i*(E(1,ii)*conj(E(2,ii)) - E(2,ii)*conj(E(1,ii)));
    
end



end

