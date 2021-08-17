%% Example 2.  IPP
% Generate data from 2-d space, let lambda(x,y)=exp(-(x-9).^2./10-(y-9).^2./10);
clear;
%addpath('C:\Users\hongxiao\Documents\teaching\stat6474_2017\lectures\Lecture_22_23_Gaussian_process_regression\Lec2_GP_application\new_organized_code')
s = RandStream('mt19937ar','Seed',152);
RandStream.setGlobalStream(s);

ngrid=[30,20]; % num of grid points on each dimension.
maxTx=500;
maxTy=500;
[t1,t2] = meshgrid(linspace(0,maxTx,ngrid(1)), linspace(0,maxTy,ngrid(2)));                                
%lam_max=.032; % maximum of the lam(s)
%lambda_true =lam_max*(exp(-(t1-61).^2/4000-(t2-91).^2/4000));
%lam_max=.0285; % maximum of the lam(s)
%lambda_true =lam_max*(exp(-(t1-62).^2/4000-(t2-141).^2/4000));
%lam_max=.0267; % maximum of the lam(s)
%lambda_true =lam_max*(exp(-(t1-71).^2/4000-(t2-221).^2/4000));
lam_max=.0305; % maximum of the lam(s)
lambda_true =lam_max*(exp(-(t1-205).^2/4000-(t2-103).^2/4000));

figure()
surf(t1,t2,lambda_true)
title('True lambda(s) surface');
xlabel('x')
ylabel('y')

%sig_fn = inline('.032*exp(-(x(:,1)-61).^2/4000-(x(:,2)-91).^2/4000)'); 
%sig_fn = inline('.0285*exp(-(x(:,1)-62).^2/4000-(x(:,2)-141).^2/4000)'); 
%sig_fn = inline('.0267*exp(-(x(:,1)-71).^2/4000-(x(:,2)-221).^2/4000)'); 
sig_fn = inline('.0305*exp(-(x(:,1)-205).^2/4000-(x(:,2)-103).^2/4000)'); 
T=[0,maxTx;0,maxTy];
[A,vol_T]=inhomo_poisson_process_sampler(lam_max,T,sig_fn);

figure()
imagesc(linspace(0,maxTx,ngrid(1))',linspace(0,maxTy,ngrid(2))',lambda_true)
hold on
plot(A(:,1),A(:,2),'.k','MarkerSize',20)    
hold off
title('IPP');

S=A;
csvwrite('points_4.csv',A)
%save('C:\Users\hongxiao\Documents\teaching\stat6474_2017\lectures\Lecture_22_23_Gaussian_process_regression\Lec2_GP_application\new_organized_code\simulated_testData_v1.mat','S','T','sig_fn',...
%    'lam_max','ngrid','t1','t2','lambda_true'); 
