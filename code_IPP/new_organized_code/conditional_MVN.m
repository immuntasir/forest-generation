% This function sample from a Gaussian process (at s_new) conditional on old insantiated
% values at {s}\cup{\ts}, suppose that the g_old, mu_old, K_old are
% consistent in order with the set {s}\cup{ts}.

function l_new = conditional_MVN(s_new, l_old, mu_old, K_old, mu_new, s_old, theta)

% Input:
% K_old is the prior covariance matrix for l_old. 
% mu_old is the prior mean of g_old.
% s_new: the set of event locations on which to sample l() from.
% s_old: the location of old values, which include the observed events and
%        the latent events. 
% theta: the parameters of the covariance kernel, now assumed known. 
% Output:
% l_new: a sample of l(.) on A. 
 

% Note: In this step we assume that the infinite dimension Gaussian l(.) is
% known, and we recorded the values of l(.) only at {s}\cup{ts}. Now when
% there are A more locations, we simply need to insantiate l(.) on more
% points. Note that this step is treated as prior step, so no need to
% consider the posterior distribution of l(.).

% For this code, I used the squared exponential(SE) kernel:
% K(x,y)=\theta_1^2 exp{||x-y||^2/(2theta_2^2)}, where \theta_2 determine the rate of correlation decay. 
% reference, e.g. http://www.cs.toronto.edu/~hinton/csc2515/notes/gp_slides_fall08.pdf
% reference, e.g. http://www.cs.ucl.ac.uk/staff/c.archambeau/ATML/atml_files/atml08_lect1c_gps.pdf
% reference, page 14,16,83,84 of Rasmussen&Carl Edward "Gaussian Processes
% for machine learning".
K_oldnew = se_kern_fast(theta, s_old, s_new);
K_new = se_kern_fast(theta, s_new);
n_old = size(K_old,1);

[chol_K_old,] = jitterChol(K_old);
iv_chol_K_old = solve_triu(chol_K_old,eye(n_old)); % try to use the solve_triu function of the lightspeed library.
iv_Kold = iv_chol_K_old*iv_chol_K_old';

ok = K_oldnew'*iv_Kold;
Mu_cond = mu_new + ok * (l_old - mu_old);
temp_cov = ok * K_oldnew; 
Cov_cond = K_new -(temp_cov + temp_cov')/2; % the trick (A+A')/2 simply make the matrix symmetric. 
%Cov_cond=K_new-ok*K_oldnew;

 [l_new,] = gaussianSample(1, Mu_cond', Cov_cond)';
 
% Note: from profiling test, I found that the sqrtm matrix in
% gaussianSample takes the majority of the time. So I want to replace this
% function by may be the randnorm  function of
% http://research.microsoft.com/en-us/um/people/minka/software/lightspeed/

% junk code
% [chol_Kold,pd]=chol(K_old);
% if pd==0  % not positive definite
%     inv_chol_Kold=solve_triu(chol_Kold,eye(n_old));    
%     inv_Kold=iv_chol_K_old*iv_chol_K_old';
% else
%     
%     % add a jitter value to covariance when computing inverse, refer to the function gaussianFastConditional.m of Michalis Titsias
%     %C:\Users\hongxiaozhu\Documents\Notes_and_All\cell_image\my_Notes\Elliptical_slice_sampling\ess_code\ess_code\ess_code\titsias\toolbox
%     tic;
%     %iv_chol_Kold=(K_old+jitter*K_new(1,1)*eye(n_old))\eye(n_old); 
%     iv_chol_Kold=(K_old+jitter*K_new(1,1)*eye(n_old))\eye(n_old); 
%     toc;
% end

