% This function generate samples from a homogeneous Poisson process with
% rate Omega>0, Omega is a scalar.
% Input:
% omega: the rate of the homogeneous poisson process.
% T:     the domain on which the poisson process is defined. If it is one
%        dimension, it is a vector of [0,t] the boundary of the domain. If it is 
%        two dimension, it is a 2 by 2 matrix of [t_x0,t_x1;t_y0,t_y1]. If it is
%        three dimension, it is a 3 by 2 matrix of [t_x0, t_x1; t_y0, t_y1; t_z0, t_z1] and so on. 
% Ouput:
% A:     the sample generated, each row is a sample. Columns are the
%        dimensions.
% vol_T: the volumen of T computed.

% Example:
% omega=5;
% T=[0,1;0.5,4;-1,3];

function [A,vol_T,num_A]=homo_poisson_process_sampler(omega,T)

% 1.Compute the volumn of the region. 
if isvector(T)
   vol_T=range(T);
else 
   vol_T=prod(range(T,2));
end

% 2. Generate a Poisson number with rate omega*vol(T), which gives a non-negative
% integer.
num_A = poissrnd(omega*vol_T);

% 3. Now uniformly sample num_A points from the region T.
if isvector(T)
     A = min(T) + (max(T)-min(T)).*rand(num_A,1);
else  
     A=repmat(min(T,[],2)',num_A,1)+repmat(range(T,2)',num_A,1).*rand(num_A,size(T,1));
end

 