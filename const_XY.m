%LBO晶体常数
global c;%光速
global ele_c;%ele_c=8.854187817e-12真空电容率
global P_wavelength;%泵浦光波长
global S_wavelength;%信号光波长
global I_wavelength;%闲置光波长
global P_R_index;%泵浦光折射率
global I_R_index;%闲置光折射率
global S_R_index;
global P_O_index;%寻常光折射率
global P_E_index;%异常光折射率
global S_angle;%信号光与泵浦光夹角
global P_angle;%双折射走离角
global I_angle;
global a_S;%信号光吸收系数
global a_P; %泵浦光吸收系数
global a_I; %闲置光吸收系数
%常数
%---------------------------------------------------------------------------------------------
c=2.99792e+8;
ele_c=8.8541e-12;%真空电容率
P_wavelength=532e-9;        %泵浦光波长
wvl =6.5e-9;    %光谱半极大全宽度:m
S_wavelength=1./(1/1053e-9+wvl*t/(t0*(1053e-9)^2));      %信号光波长
P_w=2*pi*c/P_wavelength;    %泵浦光频率
S_w=2*pi*c./S_wavelength;    %信号光频率
I_w=P_w-S_w;                %闲置光频率
I_wavelength=2*pi*c./I_w;    %闲置光波长 
S_R_index=sqrt(2.586179+0.013099./((S_wavelength*1e+6).^2-0.011893)-0.017968*(S_wavelength*1e+6).^2-(2.26e-4)*(S_wavelength*1e+6).^4);%信号光折射率(电场强度偏振沿y向）
I_R_index=sqrt(2.586179+0.013099./((I_wavelength*1e+6).^2-0.011893)-0.017968*(I_wavelength*1e+6).^2-(2.26e-4)*(I_wavelength*1e+6).^4);%闲置光折射率(电场强度偏振沿y向）
P_X_index=sqrt(2.454140+0.011249/((P_wavelength*1e+6)^2-0.011350)-0.014591*(P_wavelength*1e+6)^2-6.60*1e-5*(P_wavelength*1e+6)^4);
P_Y_index=sqrt(2.539070+0.012711/((P_wavelength*1e+6)^2-0.012523)-0.018540*(P_wavelength*1e+6)^2-2.00*1e-4*(P_wavelength*1e+6)^4);
P_Z_index=sqrt(2.586179+0.017968/((P_wavelength*1e+6)^2-0.011893)-0.017968*(P_wavelength*1e+6)^2-2.26*1e-4*(P_wavelength*1e+6)^4);
S_angle=0.5*pi/180;         %泵浦光与信号光波矢夹角
P_angle=-0.42*pi/180;
I_angle = -asin(S_R_index.*I_wavelength./I_R_index./S_wavelength*sin(S_angle));         %泵浦光与闲置光波矢夹角
P_R_index=(S_R_index(num/2)./S_wavelength(num/2)*cos(S_angle)+I_R_index(num/2)./I_wavelength(num/2).*cos(I_angle(num/2)))*P_wavelength;
angle=acos(sqrt((1/P_X_index^2-1/P_R_index^2)/(1/P_X_index^2-1/P_Y_index^2)));
P_R_index=(S_R_index(num/2)./S_wavelength(num/2)*cos(S_angle)+I_R_index(num/2)./I_wavelength(num/2).*cos(I_angle(num/2)))*P_wavelength;
d_eff=0.98e-12*cos(angle); %参量过程－有效非线性系数
K_con=S_w./sqrt(P_R_index.*S_R_index.*I_R_index)/c.*d_eff;
%-------------------------------------------------------------------------------------------------   
a_S=0.1; %信号光吸收系数
a_P=0.1; %泵浦光吸收系数
a_I=0.1; %闲置光吸收系数
dk=2*pi*(P_R_index/P_wavelength-S_R_index./S_wavelength*cos(S_angle)-I_R_index./I_wavelength.*cos(I_angle));
% plot(t,dk);hold on