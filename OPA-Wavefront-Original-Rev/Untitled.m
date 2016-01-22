clear all
nx=256;         %x－取样个数
ny=256;         %y－取样个数
nvar=3;        %参量方程个数
t0=3.2e-9;       %脉冲宽度：s
r0=1e-3;          %光斑直径:m
crstl_L=60e-3;   %晶体长度:m
z1=0;          %积分起点
z2=crstl_L;    %积分终点:m
s=10;               %
%初始电场强度
E_P0=8e+7;   
E_S0=8e+7;
E_I0=0;
%--------------------------------------------------------------------------
%-
dx=s*r0/nx;        %x－取样分辨率
dy=s*r0/ny;        %y－取样分辨率
x=linspace(-s/2*r0,s/2*r0,nx);%x－坐标
y=linspace(-s/2*r0,s/2*r0,ny);%y－坐标 
Exy_ph=wvf_Gn(x,y,r0,dx);

for lx=1:nx
     for ly=1:ny
%          E_P_out(:,lx,ly)=E_P0*exp(i*EXY(lx,ly)).*sqrt(exp(-x(lx)^20/(r0/2/log(2)^0.05)^20)*exp(-y(ly)^20/(r0/2/log(2)^0.05)^20).*exp(-t.^10/(t0/2/log(2)^0.1)^10))/sqrt(P_R_index);
         E_S_out(lx,ly)=E_S0*(exp(-(x(lx)^2)/(r0/2/log(2)^0.5)^2)*exp(-(y(ly)^2)/(r0/2/log(2)^0.2)^2));
         E_P_out(lx,ly)=E_P0*exp(i*10*Exy_ph(lx,ly))*(exp(-(x(lx)^2+(y(ly))^2)^5/(r0/2/log(2)^0.1)^10));%*exp(-i*P_w*(j*dt));
        %E_S_out(:,lx,ly)=E_S0*sqrt(exp(-(x(lx)^2+y(ly)^2)/(r0/2/log(2)^0.5)^2)*exp(-t.^2/(t0/2/(log(2))^0.5)^2))./sqrt(S_R_index);%*exp(-i*S_w(j)*(j*dt));
         E_I_out(:,lx,ly)=0;  
     end
end
%---------------------------------------------------------------------------
E=E_P_out;
[E,fx,fy]=xy_fft(E,x,y);
I=(abs(E)).^2;
figure(1)
plot(fx*1053e-9*1e3,I(nx/2,:)/max(max(I)))
hold on
E=E_S_out;
[E,fx,fy]=xy_fft(E,x,y);
I=(abs(E)).^2;
figure(1)
plot(fx*1053e-9*1e3,I(nx/2,:)/max(max(I)))
hold on