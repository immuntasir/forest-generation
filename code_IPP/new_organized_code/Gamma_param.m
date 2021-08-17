function [alpha1, beta1,mean1]=Gamma_param(mode1,var,plt,arg)
% according to the mode and variance, compute the two parameters of gamma.
% Note beta is the rate parameter, not the scale parameter. alpha is the
% shape parameter. 
alpha1=1+mode1^2/(2*var)+sqrt(4*var*mode1^2+mode1^4)/(2*var);
beta1=(alpha1-1)/mode1; 
mean1=beta1/alpha1;
if plt==1
    dens=NaN(length(arg),1);
    for k=1:length(arg)
        dens(k)=gampdf(arg(k),alpha1,1/beta1);
    end
    figure()
    plot(arg,dens)
    hold on
    plot([mean1,mean1],[0,max(dens)+0.1],'-k');
    plot([mode1,mode1],[0,max(dens)+0.1],'--k');
    hold off
end
end