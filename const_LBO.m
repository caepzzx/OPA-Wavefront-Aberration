%LBO晶体常数
global c;%光速
global ele_c;%ele_c=8.854187817e-12真空电容率
global P_wavelength;%泵浦光波长
global S_wavelength;%信号光波长
global I_wavelength;%闲置光波长
global P_R_index;%泵浦光折射率
global S_R_index;%信号光折射率[zzx]
global I_R_index;%闲置光折射率
global P_O_index;%寻常光折射率
global P_E_index;%异常光折射率
global S_angle;%信号光与泵浦光夹角
global P_angle;%双折射走离角
global a_S;%信号光吸收系数
global a_P; %泵浦光吸收系数
global a_I; %闲置光吸收系数
global wvl; %光谱半极大全宽度:m
global K_con0; %中心波长处耦合项参数
global K_con; %不同波长耦合项参数
%常数
%---------------------------------------------------------------------------------------------
c=2.99792e+8;
ele_c=8.8541e-12;%真空电容率
P_wavelength=532e-9;        %泵浦光波长

% %
% %计算闲频光中心波长
% P_wavelength0=P_wavelength;
S_wavelength0=1053e-9;      %[zzx]信号光中心波长
% P_w0=2*pi*c/P_wavelength0;    %泵浦光中心频率
% S_w0=2*pi*c/S_wavelength0;    %信号光中心频率
% I_w0=P_w0-S_w0;              %闲置光中心频率
% I_wavelength0=2*pi*c/I_w0;   %闲置光中心波长 

% 
% lambda0
% deltav=0.5/duration;
% deltalambdainv=deltav/3e8;
% wvl=lambda0^2*deltalambdainv;
wvl =6.5e-9;    %光谱半极大全宽度:m
R_wvl=3;    %控制考虑的波长范围
S_wavelength=linspace(S_wavelength0-R_wvl*wvl/2,S_wavelength0+R_wvl*wvl/2,nwav);      %考虑的信号光波长
P_w=2*pi*c/P_wavelength;    %泵浦光频率
S_w=2*pi*c./S_wavelength;    %信号光频率
I_w=P_w-S_w;                %闲置光频率
I_wavelength=2*pi*c./I_w;    %闲置光波长 
S_R_index=sqrt(2.586179+0.013099./((S_wavelength*1e+6).^2-0.011893)-0.017968*(S_wavelength*1e+6).^2-(2.26e-4)*(S_wavelength*1e+6).^4);%信号光折射率(电场强度偏振沿y向）
I_R_index=sqrt(2.586179+0.013099./((I_wavelength*1e+6).^2-0.011893)-0.017968*(I_wavelength*1e+6).^2-(2.26e-4)*(I_wavelength*1e+6).^4);%闲置光折射率(电场强度偏振沿y向）
P_X_index=sqrt(2.454140+0.011249/((P_wavelength*1e+6)^2-0.011350)-0.014591*(P_wavelength*1e+6)^2-6.60*1e-5*(P_wavelength*1e+6)^4);
P_Y_index=sqrt(2.539070+0.012711/((P_wavelength*1e+6)^2-0.012523)-0.018540*(P_wavelength*1e+6)^2-2.00*1e-4*(P_wavelength*1e+6)^4);
P_Z_index=sqrt(2.586179+0.017968/((P_wavelength*1e+6)^2-0.011893)-0.017968*(P_wavelength*1e+6)^2-2.26*1e-4*(P_wavelength*1e+6)^4);
S_angle=-0.5*pi/180;         %泵浦光与信号光波矢夹角
P_angle=0.42*pi/180;          %[zzx]信号光与闲置光波矢夹角
I_angle = -asin(S_R_index.*I_wavelength./I_R_index./S_wavelength*sin(S_angle));         %泵浦光与闲置光波矢夹角
P_R_index=(S_R_index(num/2)./S_wavelength(num/2)*cos(S_angle)+I_R_index(num/2)./I_wavelength(num/2).*cos(I_angle(num/2)))*P_wavelength;
angle=acos(sqrt((1/P_X_index^2-1/P_R_index^2)/(1/P_X_index^2-1/P_Y_index^2)));
d_eff=0.98e-12*cos(angle); %参量过程－有效非线性系数

%------计算耦合项在中心频率处的参数--------%
% S_R_index0=sqrt(2.586179+0.013099./((S_wavelength0*1e+6).^2-0.011893)-0.017968*(S_wavelength0*1e+6).^2-(2.26e-4)*(S_wavelength0*1e+6).^4);%信号光在中心频率处折射率(电场强度偏振沿y向）
% I_R_index0=sqrt(2.586179+0.013099./((I_wavelength0*1e+6).^2-0.011893)-0.017968*(I_wavelength0*1e+6).^2-(2.26e-4)*(I_wavelength0*1e+6).^4);%闲置光在中心频率处折射率(电场强度偏振沿y向）
% S_angle0=S_angle;
% I_angle0 = -asin(S_R_index0*I_wavelength0/I_R_index0/S_wavelength0*sin(S_angle0)); %中心频率处泵浦光与闲频光波矢夹角
% P_R_index0=(S_R_index0./S_wavelength0*cos(S_angle)+I_R_index0./I_wavelength0.*cos(I_angle0))*P_wavelength0;%泵浦光在中心频率处折射率
% 
% 
% K_con_S0=S_w0*d_eff/(c*S_R_index0.*cos(S_angle0));
% K_con_I0=I_w0*d_eff/(c*I_R_index0.*cos(I_angle0));
% K_con_P0=P_w0*d_eff/(c*P_R_index0);
% K_con0=[K_con_S0;K_con_I0;K_con_P0];

%-------------------------------------------------%

%参量过程中三波的耦合项中的参数[zzx]
K_con_S=S_w*d_eff./(c*S_R_index.*cos(S_angle));
K_con_I=I_w*d_eff./(c*I_R_index.*cos(I_angle));
K_con_P=P_w*d_eff./(c*P_R_index);
K_con={K_con_S;K_con_I;K_con_P};
%-------------------------------------------------------------------------------------------------   
a_S=0.1; %信号光吸收系数
a_P=0.1; %泵浦光吸收系数
a_I=0.1; %闲置光吸收系数
dk=2*pi*(P_R_index/P_wavelength-S_R_index./S_wavelength*cos(S_angle)-I_R_index./I_wavelength.*cos(I_angle));
% plot(t,dk);hold on