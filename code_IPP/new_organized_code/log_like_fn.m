% This function compute the loglikelihood function 
% f(G,G_til|\lambda(t)). Which is in the form of 
% \Omega^|E|exp^{-\Omega vol(T)} \prod_{i=1}^{|G|} \sigma(l)
% \prod_{j=1}^{|G_til|}(1-\sigma(l)). 
% Refer to (6) of my note, or (4) of Adams, Murray and MacKay (2011).

function logL=log_like_fn(l_all,N,Omega,vol_of_T)

M=length(l_all)-N;
l_s=l_all(1:N);
l_ts=l_all(N+1:end);
logL=(N+M)*log(Omega)-Omega*vol_of_T-sum(log(1+exp(-l_s)))-sum(log(1+exp(l_ts)));


