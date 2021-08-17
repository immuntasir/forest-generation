% This function generates samples from inhomogeneous poisson process
% (signoidal Gaussian Cox process) with intensity described in sig_fn 
% Reference: Algorithm 1 of Adams, Murray and MacKay.
% omega: the upper bound of the intensity. The lambda^* in the paper.
% T:     the domain on which the poisson process is defined. If it is one
%        dimension, it is a vector of [t1,t2] the boundary of the domain. If it is 
%        two dimension, it is a 2 by 2 matrix of [t_x0,t_x1;t_y0,t_y1]. If it is
%        three dimension, it is a 3 by 2 matrix of [t_x0, t_x1; t_y0, t_y1; t_z0, t_z1] and so on. 


function [A,vol_T]=inhomo_poisson_process_sampler(omega,T,sig_fn)

% 1.Compute the volumn of the region. 
if isvector(T)
   vol_T=range(T);
else 
   vol_T=prod(range(T'));
end

% step 1: generate a poisson number
num_A = poissrnd(omega*vol_T);

% step 2: conditional on the number ,generate uniform event locations.
if isvector(T)
   all_s=min(T)+(max(T)-min(T))*rand(num_A,1); 
else
   all_s=repmat(min(T,[],2)',num_A,1) + repmat(range(T,2)',num_A,1).*rand(num_A,size(T,1));
end

% step 3: generate GP at the locations.
%K_s_mat=se_kernel(all_s, all_s, theta);
%l_s=mvnrnd(zeros(size(all_s,1),1)',K_s_mat)';
sig_s=sig_fn(all_s);

% step 4
sellist=rand(num_A,1)<sig_s;
if isvector(T)
   A=all_s(sellist);
else
   A=all_s(sellist,:);
end
