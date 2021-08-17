%% Example 1.  HPP
Vol = 1;
lam = 30;
N = lam*Vol;
X=[rand(N), rand(N)];
figure()
plot(X(:,1),X(:,2),'o');
title = 'HPP';


%% Example 2.  IPP
% Generate data from 2-d space, let lambda(x,y)=exp(-(x-9).^2./10-(y-9).^2./10);
clear;
%addpath('C:\Users\hongxiao\Documents\teaching\stat6474_2017\lectures\Lecture_22_23_Gaussian_process_regression\Lec2_GP_application\new_organized_code')
s = RandStream('mt19937ar','Seed',152);
RandStream.setGlobalStream(s);

ngrid=[30,20]; % num of grid points on each dimension.
maxTx=11;
maxTy=11;
lam_max=1; % maximum of the lam(s)
[t1,t2] = meshgrid(linspace(0,maxTx,ngrid(1)), linspace(0,maxTy,ngrid(2)));                                
lambda_true =lam_max*(exp(-(t1-9).^2/10-(t2-1.5).^2/10) + exp(-(t1-3).^2/10-(t2-2).^2/10) + exp(-(t1-2.5).^2/10-(t2-8.5).^2/10) + exp(-(t1-8.5).^2/10-(t2-8.5).^2/10));

figure()
surf(t1,t2,lambda_true)
title('True lambda(s) surface');
xlabel('x')
ylabel('y')

sig_fn = inline('exp(-(x(:,1)-9).^2/10-(x(:,2)-1.5).^2/10) + exp(-(x(:,1)-3).^2/10-(x(:,2)-2).^2/10) + exp(-(x(:,1)-2.5).^2/10-(x(:,2)-8.5).^2/10) + exp(-(x(:,1)-8.5).^2/10-(x(:,2)-8.5).^2/10)'); 
T=[0,maxTx;0,maxTy];
[A,vol_T]=inhomo_poisson_process_sampler(lam_max,T,sig_fn);

figure()
imagesc(linspace(0,maxTx,ngrid(1))',linspace(0,maxTy,ngrid(2))',lambda_true)
hold on
plot(A(:,1),A(:,2),'.k','MarkerSize',20)    
hold off
title('IPP');

S=A;
csvwrite('test2.csv',A)
%save('C:\Users\hongxiao\Documents\teaching\stat6474_2017\lectures\Lecture_22_23_Gaussian_process_regression\Lec2_GP_application\new_organized_code\simulated_testData_v1.mat','S','T','sig_fn',...
%    'lam_max','ngrid','t1','t2','lambda_true'); 

%% Example 3.  Sample S, S_tilde, g from the Blocked Gibbs Sampler.
clear;
route = 'C:\Users\hongxiao\Documents\teaching\stat6474_2017\lectures\Lecture_22_23_Gaussian_process_regression\Lec2_GP_application\new_organized_code';
load([route, '\simulated_testData_v1.mat']);

% % rescale the domain T.
vol_T_test=prod(range(T,2));
T0=T;
scale=(vol_T_test/size(S,1))^(-1/size(T,1)); % the main idea is to let mu(T)=N, so lam^*=O(N). We therefore can set initial values lam_star0=1, and set Gamma(\alpha,\beta) priors to have mode at 1.
T=T0*scale; % T is rescaled so that mu(T)=N;
vol_T=prod(range(T,2));
%  
s=S*scale;
% clear S;
% % Generate initial S_tilde, and g_{S\cup S_tilde}
n_ts0=size(s,1);
ts0=repmat(min(T,[],2)',n_ts0,1)+repmat(range(T,2)',n_ts0,1).*rand(n_ts0,size(T,1));
s0_all=[s;ts0];
theta0=[5,1.5]';
lgtheta_param=[3,1];
K_mat0 = se_kern_fast(theta0,s0_all);
l0 = gaussianSample(1, zeros(1,size(s0_all,1)), K_mat0)';


nMc=8000; 
%profile on
fixlam=1;
lam_star0 =5; 
pred=0; 
grid=[];
gamma_prior =[];
tt=tic;
[Mlist,num_A_list, ts_list,l_list,loglik_list,lam_star_list,theta_list,lamb_pred]=...
    Block_gibbs_sampler_for_SGCP(s,ts0,l0,lam_star0,T,theta0,nMc,gamma_prior,lgtheta_param,pred,grid,fixlam);
time=toc(tt);

%save('C:\Users\hongxiao\Documents\teaching\stat6474_2017\lectures\Lecture_22_23_Gaussian_process_regression\Lec2_GP_application\simulated_test_Res.mat');

%% Result plot.
clear;
load('C:\Users\hongxiao\Documents\teaching\stat6474_2017\lectures\Lecture_22_23_Gaussian_process_regression\Lec2_GP_application\simulated_test_Res6.mat');
figure()
hist(Mlist)
title('Histogram of M')

figure()
imagesc(linspace(0,max(T(1,:)),ngrid(1))',linspace(0,max(T(2,:)),ngrid(2))',lambda_true*scale)
hold on
plot(s(:,1),s(:,2),'ok','LineWidth',2)    
hold off
title('True \lambda(s) and the sample')

% compute the histogram of the thinned event.
N=size(s,1);
xgrid=linspace(0,max(T(1,:)),60);
ygrid=linspace(0,max(T(2,:)),30);
Hist_count=zeros(length(ygrid)-1,length(xgrid)-1);
burnin=3000;
for i=burnin+1:nMc
    tsi=ts_list{i};
    for j=1:size(tsi,1)        
        id=find(abs(xgrid-tsi(j,1))==min(abs(xgrid-tsi(j,1))));
        if xgrid(id)>tsi(j,1)
            id_col=id-1;
        else
            id_col=id;
        end
        
        id1=find(abs(ygrid-tsi(j,2))==min(abs(ygrid-tsi(j,2))));
        if ygrid(id1)>tsi(j,2)
            id_row=id1-1;
        else
            id_row=id1;
        end
        Hist_count(id_row,id_col)=Hist_count(id_row,id_col)+1;        
    end
end

norm_hist=Hist_count/sum(sum(Hist_count));
figure()
imagesc(xgrid,ygrid,norm_hist)
hold on
plot(s(:,1),s(:,2),'ok','LineWidth',2);
hold off
title('Normed histogram of thinned events');

figure()
imagesc(xgrid,ygrid,-norm_hist)
hold on
plot(s(:,1),s(:,2),'ok','LineWidth',2);
hold off
title('Negative of the Normed histogram of thinned events');



figure()
subplot(2,2,1)
hist(Mlist(burnin+1:end))
title('Histogram of the M');
subplot(2,2,2)
hist(lam_star_list(burnin+1:end))
title('Histogram of the \lambda^*')
subplot(2,2,3)
hist(theta_list(burnin+1:end,1));
title('Histogram of the \theta_1');
subplot(2,2,4)
hist(theta_list(burnin+1:end,2));
title('Histogram of the \theta_2');

