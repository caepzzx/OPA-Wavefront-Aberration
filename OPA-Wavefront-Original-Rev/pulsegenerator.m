function E=pulsegenerator(x,y,t,t0,d0,mxy,mt)
%Pulse Generation
% mxy    高斯函数阶数-空间
% mt     高斯函数阶数-时间
% x,y,t  输入空间时间离散坐标
% t0     脉宽
% d0     光束直径
% 编写于2008/12/5。
num=size(t,2);
nx=size(x,2);
ny=size(y,2);
E=zeros(num,nx,ny);
[X,Y]=meshgrid(x,y);
for lt=1:num   
    E(lt,:,:)=sqrt(exp(-(X.^2+Y.^2).^mxy/(d0/(2*log(2)^(1/2/mxy)))^(2*mxy))*exp(-t(lt)^(2*mt)/(t0/2/log(2)^(1/2/mt))^(2*mt)));
end