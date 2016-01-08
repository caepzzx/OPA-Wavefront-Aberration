%% GVM in LBO crystal
num=50;                      %时间取样个数
L=60e-3;                     %Length of LBO
t0=40e-12;                     %脉冲宽度：ns
s=2;                         %窗口宽度因子
t=linspace(-s*t0,s*t0,num);  %t－坐标
c=2.99792e+8;
ele_c=8.8541e-12;           %真空电容率
                 
S_lambda0=800e-9; %Lambda0 of signal
lambda0=S_lambda0;%central wavelength
duration=12e-15;%transform bandwidth in various units
deltav=0.5/duration;
deltalambdainv=deltav/3e8;
wvl=lambda0^2*deltalambdainv; %光谱半极大全宽度:nm


P_lambda0=532e-9;%Lambda0 of pump
S_wavelength=linspace(S_lambda0-wvl*s^2,S_lambda0+wvl*s^2,1000);
P_wavelength=P_lambda0; 

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
P_angle=0.42*pi/180;
I_angle = -asin(S_R_index.*I_wavelength./I_R_index./S_wavelength*sin(S_angle));         %泵浦光与闲置光波矢夹角
P_R_index=(S_R_index(num/2)./S_wavelength(num/2)*cos(S_angle)+I_R_index(num/2)./I_wavelength(num/2).*cos(I_angle(num/2)))*P_wavelength;
angle=acos(sqrt((1/P_X_index^2-1/P_R_index^2)/(1/P_X_index^2-1/P_Y_index^2)));
P_R_index=(S_R_index(num/2)./S_wavelength(num/2)*cos(S_angle)+I_R_index(num/2)./I_wavelength(num/2).*cos(I_angle(num/2)))*P_wavelength;

diff_S_R=diff(S_R_index);
GVM=abs(S_R_index/c-S_wavelength/c*diff_S_R(round(numel(diff_S_R)/2))-P_R_index/c);
GVM=GVM(round(numel(GVM)/2));
DeltaTao=GVM*L;
Ratio_Tao=DeltaTao/t0