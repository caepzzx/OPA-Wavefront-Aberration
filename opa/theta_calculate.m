function theta=theta_calculate(alpha,lambdas,lambdap,material,var);
%% 计算相位匹配角度
%1.输入：
%      material：材料名称，字符串
%      alpha：    非共线角，单位：角度
%      lambdas: 信号光中心波长，单位：米
%      lambdap:泵浦光中心波长，单位：米
%      var:     预留接口，DKDP时是掺杂浓度，百分数
%2.输出：
%      theta：相位匹配角
%2012-11-01 郭晓杨
%2013-01-09 修改
%2013-01-16 修改
c=3e8;
alpha=alpha*pi/180;
lambdai=1./(1./lambdap-1./lambdas);



if nargin==4
    no_s=sellmeier(lambdas,material,'o');
ne_s=sellmeier(lambdas,material,'e');

no_i=sellmeier(lambdai,material,'o');
ne_i=sellmeier(lambdai,material,'e');

no_p=sellmeier(lambdap,material,'o');
ne_p=sellmeier(lambdap,material,'e');

elseif nargin==5;
    
    no_s=sellmeier(lambdas,material,'o',var);
ne_s=sellmeier(lambdas,material,'e',var);

no_i=sellmeier(lambdai,material,'o',var);
ne_i=sellmeier(lambdai,material,'e',var);

no_p=sellmeier(lambdap,material,'o',var);
ne_p=sellmeier(lambdap,material,'e',var);

else
    errordlg('输入参数个数有误','theta_calculate')
end


ks=2*pi*no_s./lambdas;
ki=2*pi*no_i./lambdai;
gamma=asin(ks./ki.*sin(alpha));
kp=ks.*cos(alpha)+ki.*cos(gamma);
sintheta=ne_p./kp.*sqrt(((2*pi./lambdap).^2.*no_p.^2-kp.^2)./(no_p.^2-ne_p.^2));
theta=asin(sintheta);
theta=theta*180/pi;




