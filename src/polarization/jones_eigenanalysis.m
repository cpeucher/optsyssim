function [lambda,V1,V2,dtau] = jones_eigenanalysis(T,df)
% Jones eigenanalysis for differential group delay characterisation
%
% Parameters
% ==========
% T : complex 2x2xN array
%     Jones matrices at a frequency grid with separation df
%     T(:,:,ii) is the Jones matrix at frequency indexed by ii
%
% df : real scalar
%      frequency step at which the input Jones matrices are provided, in Hz
%
% Returns
% =======
% lambda : real 2xN array
%          Eigenvalues of T(w+dw)T^(-1)(w)
% 
% V1 : complex 2xN array
%      Eigenvectors corresponding to eigenvalues lambda(1,:) 
%
% V2 : complex 2xN array
%      Eigenvectors corresponding to eigenvalues lambda(2,:) 
%
% dtau : real 1xN array
%        Differential group delay, in s 
%        If omega is the frequency grid, dtau is evaluated at
%        omega(1:N-1).


[~,~,nsamples] = size(T);

V1 = zeros(2,nsamples - 1);
V2 = zeros(2,nsamples - 1);
lambda = zeros(2,nsamples -1);
dtau = zeros(1,nsamples - 1);

for ifreq = 1:nsamples - 1
    
    NN = T(:,:,ifreq + 1)*inv(T(:,:,ifreq));
    
    [VV,BB] = eig(NN);
    
    lambda(:,ifreq) = diag(BB);
    
    V1(:,ifreq) = VV(:,1);
    V2(:,ifreq) = VV(:,2);
    
    dtau(ifreq) = abs(angle(lambda(1,ifreq)/lambda(2,ifreq)))/df/2/pi;
    
end

end