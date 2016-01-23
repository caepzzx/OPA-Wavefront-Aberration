function E=pulsegenerator(x,y,t,t0,d0,mxy,mt,chirp)
%Pulse Generation
% mxy    高斯函数阶数-空间
% mt     高斯函数阶数-时间
% x,y,t  输入空间时间离散坐标
% t0     脉宽
% d0     光束直径
% 编写于2008/12/5。

c=2.99792e+8; %光速
ele_c=8.8541e-12;%真空电容率
num=size(t,2);
nx=size(x,2);
ny=size(y,2);
E=zeros(num,nx,ny);
[X,Y]=meshgrid(x,y);
T0=t0/(2*log(2)^(1/2/mxy));
D0=d0/(2*log(2)^(1/2/mxy));
dt=t(2)-t(1);
for lt=1:num   
    E(lt,:,:)=exp(-(X.^2+Y.^2).^mxy/D0^(2*mxy))*exp(-0.5*(1+1i*chirp)*(t(lt)/T0)^(2*mt));
end