function yout=rk4(y,z,h,P_w,S_w,I_w,K_con,dk)
%y－初值
%n－方程个数
%z—传播距离
%h—积分步长
%经典四阶Runge-Kutta法求微分方程
%dy/dx=f(x,y) ,x belong to [a,b]
%y(x0)=y0
%的数值解
%---------------
%Yk+1=Yk+h/6*(K1+2K2+2K3+K4)
%K1=f(Xk,Yk)
%K2=f(Xk+h/2,Yk+h/2*K1)
%K3=f(Xk+h/2,Yk+h/2*K2)
%K4=f(Xk+h,Yk+h*K3)
%RK45 此程序主要计算入射光,经过一个切片由于不同性质的光之间的耦合,而产生的新的光场分布
%----------------

n=size(y);
dydz=zeros(n);
yt=zeros(n);
dyt=yt;
dym=yt;
dy4=yt;
yout=yt;
hh=h*0.5;
h6=h/6.0;
zh=z+hh;
dydz=derivs(z,y,P_w,S_w,I_w,K_con,dk);%对应k1
yt=y+hh.*dydz;
dyt=derivs(zh,yt,P_w,S_w,I_w,K_con,dk);%对应k2
yt=y+hh.*dyt;
dym=derivs(zh,yt,P_w,S_w,I_w,K_con,dk);%对应k3
yt=y+h.*dym;
dy4=derivs(z+h,yt,P_w,S_w,I_w,K_con,dk);%对应k4
yout=y+h6.*(dydz+2.0.*dyt+2.0.*dym+dy4);%输出积分值