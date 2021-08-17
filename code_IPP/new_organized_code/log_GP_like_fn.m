% This function compute the log GP likelihood function 
% l_all \sim GP(l_all|mu_l, K_h).
% where h=[h1,h2] is the log of the two parameters of the covariance kernel.
% h1=log(theta(1)), h2=log(theta(2)). Note: theta is defined as
% K(x1,x2)=theta(1)^2 exp{-||x1-x2||^2/theta(2)^2/2}
% Input:
% h: log of theta. 
% mu_l: the mean of the l_{s}U{ts} in the GP prior. 
% l_all: current values of the l_{s}U{ts}.
% kern_fun: the handle of the kernel function. 

function [log_GP,chol_Kh]=log_GP_like_fn(h, mu_l, l_all, kern_fun, s_all)

ns=size(s_all,1);

Kh=kern_fun(exp(h),s_all);

[chol_Kh, er] = jitterChol(Kh);
if er>0
    warning('jitterChol fails, used sqrtm');
    chol_Kh=xsqrtm(Kh); 
    inv_L_t=chol_Kh\eye(ns);
    
else
    %inv_L_t=chol_Kh'\eye(ns);    
    inv_L_t=solve_triu(chol_Kh,eye(ns))';
end

tempT=inv_L_t*(l_all-mu_l);
log_GP=-ns*log(2*pi)/2-real(sum(log(diag(chol_Kh))))-sum(tempT.^2)/2;
