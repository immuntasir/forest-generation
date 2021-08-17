% This is the function to perform Bayesian Inference in Sigmoidal Gaussian 
% Cox Process, using the Algorithm 1 of Rao and Teh 2011, "Gaussian process 
% modulated renewal processes".

% Input:
% (1)s:     {s_i, i=1,\dots, n} are the observed locations of an inhomogenuous
%            Poisson process. These are data. n by dim matrix.  
% (2)ts0:    Initial values of {ts_k, k=1,\dots, M} are the latent thinned
%            events. Note that M is also a latent variable. 
% (3)l0:     Initial values of l({s_i}\cup {ts_k}), the Gaussian processes
%            instantiated at the observed events {s_i} and latent events.
% (4)lam_star: the upperbound of lambda(t). Assumed known here. Later can
%            add Inverse-Gamma prior to it.
% (5)T:      the domain on which the poisson process is defined. If it is one
%            dimension, it is a vector of [0,t] the boundary of the domain. If it is 
%            two dimension, it is a 2 by 2 matrix of [t_x0,t_x1;t_y0,t_y1]. If it is
%            three dimension, it is a 3 by 2 matrix of [t_x0, t_x1; t_y0, t_y1; t_z0, t_z1] and so on.
% (6)theta:  the parameter for the covariance kernel in GP prior. 
% (7)nMc:    Total number of iterations.
% (8)gamma_prior: a 2 by 1 vector, (\alpha,\beta) is the prior parameter of
%            lamb_star. Note that the beta is the rate parameter, not the
%            scale parameter. See wiki's definition of gamma distribution. 
% lgtheta_param: the standard deviations on the diagonal for prior variance
%             of log(theta), assuming logtheta ~ N(0,
%             diag(lgtheta_param(1)^2,lgtheta_param(2)^2));
% (8)pred:   1/0, whether to predict the lambda values on a prespecified
%            grid. 
% (9)grid:   If pred=1, a grid on T over which the values of lambda are predicted on in
%            each iteration. grid should be in the same format as s.
% fixlam: 1/0 whether fix lambda.

% Output:
% Mlist: the list of number of thinned events. 
% num_A_list: the total number of latent events generated.
% ts_list: the list of tilde{s} values across MCMC.
% l_list: the l_{s}U{ts} values.
% loglike_list: the log likelihood values (excluding the GP prior and prior
%          for h).
% lam_star_list: the list of lam_star values.
% theta_list: the list of theta values.
% lamb_pred: the predictive lambda values on grid. 


function [Mlist,num_A_list, ts_list,l_list,loglik_list,lam_star_list,theta_list,lamb_pred]=Block_gibbs_sampler_for_SGCP(s,ts0,l0,lam_star0,T,theta0,nMc,gamma_prior,lgtheta_param,pred,grid,fixlam)


if pred == 0 && nargin < 9
    grid = [];
end

l_all = l0;
lam_star = lam_star0;
ts = ts0;
theta=theta0;
N = size(s,1);
K_ss = se_kern_fast(theta, s);
K_st = se_kern_fast(theta, s, ts);
K_tt = se_kern_fast(theta,ts);
K_all = [K_ss,K_st;K_st',K_tt];    
s_all = [s;ts];

num_A_list = NaN(nMc, 1);
Mlist = NaN(nMc, 1);
ts_list = cell(nMc, 1);
l_list = cell(nMc, 1);
loglik_list = NaN(nMc, 1);
if fixlam == 0
  lam_star_list = NaN(nMc, 1);
end
theta_list = NaN(nMc, length(theta));
if pred == 1
    lamb_pred = NaN(nMc,size(grid, 1));
else
    lamb_pred = [];
end

log_GP_like_h = @log_GP_like_fn;
kern_fun_h = @se_kern_fast;
    

for i=1:nMc
    
    %sprintf(['Iteration ',int2str(i)])
    %sprintf('Iter=%d, lam_star = %10.2f',i, lam_star)
                
    %% Step 1. Sample A from a homogeneous poisson process with rate lams
    % (lambda^* in derivation). 
    [A, vol_T, numA] = homo_poisson_process_sampler(lam_star, T);     
    num_A_list(i) = numA;
        
    %% Step 2. Sample g_A|g_old from the conditional multivariate normal.     
    mu_new = zeros(numA, 1);
    mu_old = zeros(size(s_all,1), 1);       
    %profile on;
    %tic
    l_A = conditional_MVN(A, l_all, mu_old, K_all, mu_new, s_all, theta); % note: it is really slow, especially when A is large. 
    %tim=toc;  
    %profile viewer;
    %% Step 3. Thinning, keep the elements of A in interval [G_{i-1}, G_i] with
    % probability [1-\sigma(l)]. Note, for Poisson process case, since [1-sigma(l)] does not depend on G_i or
    % G_{i-1}, we can do this all together using matrix operation.
    thres = 1 - 1./(1 + exp(-l_A));  % logistic transformation of l_A.
    selA = (rand(numA,1) < thres);
    ts = A(selA, :);  % the new tilde{s}.
    ts_list{i} = ts;
    M = size(ts,1);
    Mlist(i) = M;   
    s_all = [s;ts];
    K_st = se_kern_fast(theta, s, ts);
    K_tt = se_kern_fast(theta, ts);
    K_all = [K_ss, K_st; K_st', K_tt]; % note that K_ss also need to be updated when theta is updated.       
    l_all = [l_all(1:N); l_A(selA)];
     
    sprintf('Iter=%d, num_A=%d, M=%d', i, numA, M)
    
    %% Step 4. Updating l_{g_{{s_i}\cup{ts_m}}. 
    cur_log_like = log_like_fn(l_all, N, lam_star, vol_T);
    log_like_fn_h = @log_like_fn;
    [nu_sample, chol_K_all, err] = gaussianSample(1, zeros(1,length(l_all)), K_all);
    [l_all, cur_log_like] = elliptical_slice(l_all, nu_sample', log_like_fn_h, cur_log_like, 0, N, lam_star, vol_T); % l_all is updated. 
    l_list{i} = l_all;
    loglik_list(i) = cur_log_like;
    if pred == 1        
        mu_old = zeros(size(s_all, 1), 1);
        l_pred = conditional_MVN(grid, l_all, mu_old, K_all, zeros(size(grid, 1), 1), s_all, theta);
        lamb_pred(i,:) = lam_star./(1 + exp(-l_pred))';
    end
    
    
    %% Step 5. Update lambda^* from the Gamma distribution. Based on my test, the mixing is bad when this step is added.
    if fixlam==0
       lam_star = gamrnd(M + N + gamma_prior(1), 1 /(vol_T + gamma_prior(2)));  
       lam_star_list(i) = lam_star;
    end
    
    %% Step 6. Update theta parameters.
    % compute the colesky decomponsition of K_all.    
    if err>0              
        temp_quad = (chol_K_all\eye(length(l_all)))*l_all;
    else
        temp_quad = solve_triu(chol_K_all,eye(length(l_all)))'*l_all;
    end
    cur_log_like_GP = -(N+M)*log(2*pi)/2-real(sum(log(diag(chol_K_all))))-max(sum(temp_quad.^2),0)/2; % in case sum(temp_quard.^2) is complex, we take 0. 
    [nu_theta,] = gaussianSample(1, zeros(1,length(theta)), diag(lgtheta_param.^2));
    [lg_theta,] = elliptical_slice(log(theta), nu_theta', log_GP_like_h, cur_log_like_GP, 0, zeros(length(l_all), 1), l_all, kern_fun_h, s_all); % h is updated.     
    theta = exp(lg_theta); % theta is updated here
    K_ss = se_kern_fast(theta, s); % K_ss need to be updated here whenever theta is updated.
    K_st = se_kern_fast(theta, s, ts);
    K_tt = se_kern_fast(theta,ts);
    K_all = [K_ss,K_st;K_st',K_tt];    
    theta_list(i,:) = theta';
            
end



