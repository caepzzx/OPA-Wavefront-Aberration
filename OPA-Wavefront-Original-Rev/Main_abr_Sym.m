% function main
tic
clear all;
%步长初始值设置
%----------------------------------------------  
num=50;                      %时间取样个数
nstep=50;                    %z－积分步长个数
nx=128;                      %x－取样个数
ny=128;                      %y－取样个数
nvar=3;                      %参量方程个数
t0=3e-9;                     %脉冲宽度：ns
wvl =8e-9;                   %光谱半极大全宽度:nm
d0=2.0e-3;                   %光斑直径:m
crstl_L=59.5e-3;             %晶体长度:m
z1=0;                        %积分起点
z2=crstl_L;                  %积分终点:m
s=2;                         %窗口宽度因子
%初始电场强度
E_P0=8.9e+7;   
E_S0=5.0e+3;
E_I0=0;
%---------------------------------------------------------------------------
dt=s*t0/num;                %时间取样分辨率
dx=s*d0/nx;                 %x－取样分辨率
dy=s*d0/ny;                 %y－取样分辨率
x=linspace(-s*d0,s*d0,nx);  %x－坐标
y=linspace(-s*d0,s*d0,ny);  %y－坐标 
t=linspace(-s*t0,s*t0,num); %t－坐标
%全局变量
%-----------------
const_LBO;         
%----------------
zz=zeros(nstep+1,1);
E_P_out=zeros(num,nx,ny);   %泵浦光电场强度
E_S_out=zeros(num,nx,ny);   %信号光电场强度
E_I_out=zeros(num,nx,ny);   %闲置光电场强度
E_P_w=zeros(num,nx,ny);     %泵浦光电场频域表示
E_S_w=zeros(num,nx,ny);     %信号光电场频域表示
E_I_w=zeros(num,nx,ny);     %闲置光电场频域表示
E_S_ph=zeros(nstep+1,nx,ny);    %信号光波前在参量作用过程中的变化
E_I_ph=zeros(nstep+1,nx,ny);    %闲置光波前在参量作用过程中的变化
E_P_ph=zeros(nstep+1,nx,ny);    %泵浦光波前在参量作用过程中的变化
E_S_am=zeros(nstep+1,nx,ny);    %信号光电场强度在参量作用过程中的变化
E_I_am=zeros(nstep+1,nx,ny);    %闲置光电场强度在参量作用过程中的变化
E_P_am=zeros(nstep+1,nx,ny);    %泵浦光电场强度在参量作用过程中的变化
E_S_xy=zeros(nx,ny);
E_I_xy=zeros(nx,ny);
E_P_xy=zeros(nx,ny);
Bmxz=zeros(nstep+1,num);
Bmxy=zeros(nstep+1,num);
bmx=zeros(1,num);
bmy=zeros(1,num);
v=zeros(nvar,nx,ny);
sum_Es=zeros(nstep+1,1);
sum_Ep=zeros(nstep+1,1);
Is=zeros(num,nx,ny);
Ip=zeros(num,nx,ny);
P=zeros(num,1);
[X,Y]=meshgrid(x,y);
% EXY=normrnd(1,0.0625,nx,ny);% EXY=cos(15*(X+Y)/d0)*pi*2;% EXY=normrnd(1,0.0625,nx,ny);%引入随机噪声
%畸变波前产生函数
%-------------------------------------------------------------------------
% Exy_ph=wvf_Gn(x,y,d0,dx);
% save('data\ph_abr2.mat','Exy_ph');
buf=load('data\ph_abr2.mat');
Exy_ph=buf.Exy_ph;
%-------------------------------------------------------------------------
W_F=exp(-(X.^2+Y.^2).^20/(1.3*d0/2/log(2)^(0.025))^40);
%电场强度赋初值
%--------------------------------------------------------------------------
E_S_out=E_S0*pulsegenerator(x,y,t,t0,d0,1,1)/sqrt(S_R_index(num/2));
E_P_out=E_P0*pulsegenerator(x,y,t,t0,d0,5,5)/sqrt(P_R_index);%.*exp(i*0.6*Exy_ph);
%--------------------------------------------------------------------------
for k=1:nstep
    E_P_out(k,:,:)=exp(i*0.6*Exy_ph).*squeeze(E_P_out(k,:,:));
end
%画出信号光、闲置光初始时间波形
figure(1) 
subplot(2,2,1)
I=(1/2*c*ele_c).*P_R_index*E_P_out(:,nx/2,ny/2).*conj(E_P_out(:,nx/2,ny/2));
plot(t*1e9,I/max(I),'k','LineWidth',1);
xlabel('time (ns)','FontSize',16);ylabel('Normalized intensity','FontSize',16);
hold on;
I=(1/2*c*ele_c).*E_S_out(:,nx/2,ny/2).*conj(E_S_out(:,nx/2,ny/2));
plot(t*1e9,I/max(I),'r','LineWidth',1);

%泵浦光、信号光初始能量
%----------------------------------------------------
Is=(1/2*c*S_R_index(num/2)*ele_c).*E_S_out.*conj(E_S_out);
Ip=(1/2*c*P_R_index*ele_c).*E_P_out.*conj(E_P_out);
S_En=trapz(y,squeeze(trapz(x,squeeze(trapz(t,Is,1)),1)))*1000;
P_En=trapz(y,squeeze(trapz(x,squeeze(trapz(t,Ip,1)),1)))*1000;
%----------------------------------------------------

%由龙格－库塔法计算波形沿Z轴传播的变化
%--------------------------------------------------------------------------

h=(z2-z1)/nstep;
z=z1+h/2;
zz(1)=h/2;
[fx,fy]=spati_vector(nx,ny,x,y);
[FX,FY]=meshgrid(fx,fy); 
 parfor j=1:num
     %------------------------------------------------------------------
        %第二步计算参量作用过程
        v=[E_S_out(j,:,:);E_I_out(j,:,:);E_P_out(j,:,:)]; 
        v=rk4(v,0,ny,h/2,P_w,S_w(j),I_w(j),K_con(j),dk(j));
        E_S_out(j,:,:)=v(1,:,:);
        E_I_out(j,:,:)=v(2,:,:);
        E_P_out(j,:,:)=v(3,:,:);      
end
 [E_S_am(1,:,:),E_I_am(1,:,:),E_P_am(1,:,:),E_S_ph(1,:,:),E_I_ph(1,:,:),E_P_ph(1,:,:),sum_Ep(1),sum_Es(1)]=Recod_data(E_S_out,E_I_out,E_P_out,num,x,y,t,t0);

  
for k=2:nstep+1   
%     if z>3e-2
%        P_angle=-0.42*pi/180;
%        S_angle=0.5*pi/180;
%        I_angle(:)= -asin(S_R_index.*I_wavelength./I_R_index./S_wavelength*sin(S_angle));
%     end
    parfor j=1:num
        %第一步计算走离效应的影响
        %------------------------------------------------------------------
        %进行傅立叶变换
        E_S_xy=squeeze(E_S_out(j,:,:));
        E_I_xy=squeeze(E_I_out(j,:,:));
        E_P_xy=squeeze(E_P_out(j,:,:));
        
        [E_S_w(j,:,:),~,~]=xy_fft(E_S_xy,x,y);
        [E_I_w(j,:,:),~,~]=xy_fft(E_I_xy,x,y);
        [E_P_w(j,:,:),~,~]=xy_fft(E_P_xy,x,y);
    end
       
        E_S_w=ifft(E_S_w,[],1);
        E_I_w=ifft(E_I_w,[],1);
        E_P_w=ifft(E_P_w,[],1);
        
    parfor j=1:num
        E_S_xy=squeeze(E_S_w(j,:,:));
        E_I_xy=squeeze(E_I_w(j,:,:));
        E_P_xy=squeeze(E_P_w(j,:,:)); 
        E_S_xy = E_S_xy.*exp(((S_angle*i*2*pi)*FY-1/2*a_S)*h).*exp(-i*pi*S_wavelength(num/2)/S_R_index(num/2)*(FX.^2+FY.^2)*h);      
        E_I_xy = E_I_xy.*exp(((I_angle(num/2)*i*2*pi)*FY-1/2*a_I)*h).*exp(-i*pi*I_wavelength(num/2)/I_R_index(num/2)*(FX.^2+FY.^2)*h);     
        E_P_xy = E_P_xy.*exp(((P_angle*i*2*pi)*FY-1/2*a_P)*h).*exp(-i*pi*P_wavelength/P_R_index*(FX.^2+FY.^2)*h);
        E_S_w(j,:,:)=E_S_xy(:,:);
        E_I_w(j,:,:)=E_I_xy(:,:);
        E_P_w(j,:,:)=E_P_xy(:,:);
    end
    
    
     %----------------进行时间频率域反变换---------------%        
        E_S_w=fft(E_S_w,[],1);
        E_I_w=fft(E_I_w,[],1);
        E_P_w=fft(E_P_w,[],1);
      %------------进行空间傅立叶逆变换-----------------%  
     parfor j=1:num
        E_S_xy=squeeze(E_S_w(j,:,:));
        E_I_xy=squeeze(E_I_w(j,:,:));
        E_P_xy=squeeze(E_P_w(j,:,:));
        
        E_S_xy=xy_ifft(E_S_xy,x,y);
        E_I_xy=xy_ifft(E_I_xy,x,y);
        E_P_xy=xy_ifft(E_P_xy,x,y);
        
        E_S_out(j,:,:)=E_S_xy(:,:);
        E_I_out(j,:,:)=E_I_xy(:,:);
        E_P_out(j,:,:)=E_P_xy(:,:);  
     end
  
        %------------------------------------------------------------------
        %第二步计算参量作用过程
        parfor j=1:num
        v=[E_S_out(j,:,:);E_I_out(j,:,:);E_P_out(j,:,:)]; 
        v=rk4(v,z,ny,h,P_w,S_w(j),I_w(j),K_con(j),dk(j));
        E_S_out(j,:,:)=v(1,:,:);
        E_I_out(j,:,:)=v(2,:,:);
        E_P_out(j,:,:)=v(3,:,:);  
        end
%     [bmx,bmy]=Beam_Quality(x,y,fx,fy,S_wavelength(num/2),E_P_out);
%     Bmxz(k,:)=bmx;
%     Bmyz(k,:)=bmy;
   % 记录数据
    [E_S_am(k,:,:),E_I_am(k,:,:),E_P_am(k,:,:),E_S_ph(k,:,:),E_I_ph(k,:,:),E_P_ph(k,:,:),sum_Ep(k),sum_Es(k)]=Recod_data(E_S_out,E_I_out,E_P_out,num,x,y,t,t0);
   % 下一次迭代
    z=z+h;
    zz(k)=z; 
end


    parfor j=1:num
        v=[E_S_out(j,:,:);E_I_out(j,:,:);E_P_out(j,:,:)]; 
        v=rk4(v,z,ny,-h/2,P_w,S_w(j),I_w(j),K_con(j),dk(j));
        E_S_out(j,:,:)=v(1,:,:);
        E_I_out(j,:,:)=v(2,:,:);
        E_P_out(j,:,:)=v(3,:,:);      
    end

[E_S_am(nstep+1,:,:),E_I_am(nstep+1,:,:),E_P_am(nstep+1,:,:),E_S_ph(nstep+1,:,:),E_I_ph(nstep+1,:,:),E_P_ph(nstep+1,:,:),sum_Ep(nstep+1),sum_Es(nstep+1)]=Recod_data(E_S_out,E_I_out,E_P_out,num,x,y,t,t0);
z=z-h/2;
zz(nstep+1)=z; 
%--------------------------------------------------------------------------
%画出放大后光斑中心处电场强度的时间波形
figure(1)
subplot(2,2,1)
I=(abs(E_S_out(:,nx/2,ny/2))).^2*(1/2*c*ele_c).*S_R_index(:);
plot(t*1e9,I/max(I),'b','LineWidth',1);
%画出光斑中心处电场强度随Z轴的变化
subplot(2,2,2)
plot(zz(1:nstep+1)*1000,sum_Es,'r','LineWidth',1); hold on
plot(zz(1:nstep+1)*1000,sum_Ep,'b','LineWidth',1);
xlabel('LBO Crystal Length：mm');ylabel('Output Energy:mJ');
legend('Signal Wave','Pump Wave');
hold on;
for lx=1:1:nx
    for ly=1:1:ny
        I_S=(abs(E_S_out(:,lx,ly))).^2*(1/2*c*ele_c).*S_R_index(:);
        Z(lx,ly)=trapz(t,I_S);
    end
end
Z=Z/max(max(Z));
%画出光斑形状
subplot(2,2,3)
pcolor(x,y,Z);
colormap jet,shading interp;
colorbar;
title('时间积分光斑形状');
pha=phase(E_S_out(:,nx/2,ny/2))/pi;
subplot(2,2,4)
plot(t*1e9,pha,'r');
xlabel('time (ns)');
ylabel('Phase (pi)');
title('Phase accumulated in OPA');
hold on
%画出参量作用后波前变化
figure(2)
subplot(2,2,1)
title('空间强度调制');
mesh(X*1e3,Y*1e3,Z)
hold on
subplot(2,2,2)
E(:,:)=E_S_ph(nstep+1,:,:);
mesh(X*1e3,Y*1e3,E.*W_F)
subplot(2,2,3)
E(:,:)=E_I_ph(nstep+1,:,:);
mesh(X*1e3,Y*1e3,E.*W_F)
subplot(2,2,4)
E(:,:)=E_P_ph(nstep+1,:,:);
mesh(X*1e3,Y*1e3,E.*W_F)
% 保存信号场
save('data\E_S_out.mat','E_S_out');

toc    
