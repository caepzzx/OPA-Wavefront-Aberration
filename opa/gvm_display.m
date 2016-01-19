clear all
close all
clc
profile on
tic
alpha=0;
N=1000;
material='BBO';
matchtype='ooe';
lamdap=.8e-6;
lamdas=linspace(1e-6,1.6e-6,N);
theta=pma_theta(N,alpha,lamdas,lamdap,material,matchtype);
[delta_ip1,delta_sp1,alpha1]=gvm_delta(N,lamdas,lamdap,theta,material,'ooe');
[delta_ip2,delta_sp2,alpha2]=gvm_delta(N,lamdas,lamdap,theta,material,'oee');
[delta_ip3,delta_sp3,alpha3]=gvm_delta(N,lamdas,lamdap,theta,material,'eoe');
plot(lamdas*1e6,delta_ip1*1e12,'rd',lamdas*1e6,delta_sp1*1e12,'gd')
hold on
plot(lamdas*1e6,delta_ip2*1e12,'b',lamdas*1e6,delta_sp2*1e12,'c')
hold on
plot(lamdas*1e6,delta_ip3*1e12,'m',lamdas*1e6,delta_sp3*1e12,'y')

xlim([1.0,1.6])
ylim([-100,50])
legend('os+oi<-ep,\delta_p_i','os+oi<-ep,\delta_p_s','os+ei<-ep,\delta_p_i','os+ei<-ep,\delta_p_s','es+oi<-ep,\delta_p_i','os+ei<-ep,\delta_p_s',4)
title('group velocity mismaching')
xlabel('\lambda_s[um]')
ylabel('\delta_p_i[fs/mm]')
toc