function deltak=deltak_calculate(alpha,theta,lambdas,lambdap,material,var);
%typeI 类相位匹配  相位失配量计算
%1.输入：
%      alpha    ： 非共线夹角，度
%      theta    :相位匹配角，度
%      lambdas:  信号光波长，米
%      lambdap： 泵浦光波长，米
%      material: 材料名称，字符串
%      var:      预留接口，DKDP时为掺杂浓度
%2.输出：
%      deltak:   相位失配量
% 2012-01-08 郭晓杨
%2013-01-16 修改
c=3e8;
diel=8.854e-12; 

alpha=alpha*pi/180;
theta=theta*pi/180;
lambdai=1./(1./lambdap-1./lambdas);

if nargin==5
    no_s=sellmeier(lambdas,material,'o');
ne_s=sellmeier(lambdas,material,'e');

no_i=sellmeier(lambdai,material,'o');
ne_i=sellmeier(lambdai,material,'e');

no_p=sellmeier(lambdap,material,'o');
ne_p=sellmeier(lambdap,material,'e');

elseif nargin==6;
    
no_s=ellmeier(lambdas,material,'o',var);
ne_s=sellmeier(lambdas,material,'e',var);

no_i=sellmeier(lambdai,material,'o',var);
ne_i=sellmeier(lambdai,material,'e',var);

no_p=sellmeier(lambdap,material,'o',var);
ne_p=sellmeier(lambdap,material,'e',var);

else
    errordlg('输入参数个数有误','theta_calculate')
end



np=sqrt(1./(cos(theta).^2./no_p.^2+sin(theta).^2./ne_p.^2));%泵浦光折射率
kp=2*pi.*np./lambdap;

ks=2*pi*no_s./lambdas;

kie=2*pi*no_i./lambdai;%由能量守恒得到的波矢,e代表能量energy
kip=sqrt(ks.^2+kp.^2-2*ks.*kp.*cos(alpha));%由动量守恒得到的波矢,p代表动量

deltak=kip-kie;

