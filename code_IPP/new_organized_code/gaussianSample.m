function [X,L,r] = gaussianSample(N, mu, Sigma)
%Description:  It generates N samples from a d-dimensional
%              Gaussian distribution.
%
%N  : Number of samples to generate
%mu : mean, should be 1 by D vector. 
%Sigma : covariance matrix for the samples, should be positive definite
% zhu: note: X is N by d matrix, mu is 1 by d, a row.
d = size(Sigma,1);

X=randn(N,d); 

% Cholesky decomposition for a positive definite Covariance
% If the covariance is not positive definite uses the square root 
[L,r] = jitterChol(Sigma);
if r>0 
    %warning('r>0, using sqrtm in gaussianSampling.');
    L=xsqrtm(Sigma);
    X = real(X*L); %[X, arg2, condest] = xsqrtm(A) my own version, replace some code by c code.
else
    X = X*L;
end

X=X+repmat(mu,N,1);
