%计算光束质量
tic
I=E_S_am.*conj(E_S_am);
%计算质心坐标
for k=1:num
    I_nlz(k)=trapz(x,(trapz(y,squeeze(I(k,:,:)),2)));
    y_av(k)=trapz(x,trapz(y,Y.*squeeze(I(k,:,:)),2))/I_nlz(k);
end
for k=1:num
    E_fft(k,:,:)=fftshift(fft2(squeeze(E_S_am(k,:,:))));
end
I_fft=E_fft.*conj(E_fft);
%计算倾斜角
for k=1:num
    y_beta_nzl(k)=trapz(fy,squeeze(trapz(fx,squeeze(I_fft(k,:,:)),1)));
    y_beta(k)=-S_wavelength(num/2)*trapz(fy,squeeze(trapz(fx,FY.*squeeze(I_fft(k,:,:)),2)))/y_beta_nzl(k);
end
%计算Ay
for k=1:num
    y_d_Efft(k,:,:)=-i*2*pi*fftshift(fft2(squeeze(E_S_am(k,:,:)).*Y)).*conj(squeeze(E_fft(k,:,:)));
    y_d_Efft(k,:,:)=-y_d_Efft(k,:,:)+conj(y_d_Efft(k,:,:));
    y_A(k)=i*S_wavelength(num/2)/(2*pi)*trapz(fy,squeeze(trapz(fx,FY.*squeeze(y_d_Efft(k,:,:)),1)))/y_beta_nzl(k);
end
%计算空间方差
for k=1:num
    y_theta(k)=trapz(y,squeeze(trapz(x,(Y-y_av(k)).^2.*squeeze(I(k,:,:)),1)))/I_nlz(k);
end
%计算空间频率方差
for k=1:num
    fy_theta(k)=trapz(fx,(trapz(fy,(S_wavelength(num/2)*FY-y_beta(k)).^2.*squeeze(I_fft(k,:,:)),2)))/S_wavelength(num/2)^2/y_beta_nzl(k);
end
%计算Zy
Zy=(y_A+2*y_beta.*y_av)./(2*S_wavelength(num)^2*fy_theta);
y_theta_0=y_theta-Zy.^2*S_wavelength(num)^2.*fy_theta;
y_M_Squ=4*pi*y_theta_0.^0.5.*fy_theta.^0.5;
plot(zz(2:num+1)*1000,y_M_Squ)
hold on

%X向M平方因子
%计算光束质量
I=E_S_am.*conj(E_S_am);
%计算质心坐标
for k=1:num
    I_nlz(k)=trapz(x,(trapz(y,squeeze(I(k,:,:)),2)));
    x_av(k)=trapz(y,trapz(x,X.*squeeze(I(k,:,:)),1))/I_nlz(k);
end
for k=1:num
    E_fft(k,:,:)=fftshift(fft2(squeeze(E_S_am(k,:,:))));
end
I_fft=E_fft.*conj(E_fft);
%计算倾斜角
for k=1:num
    x_beta_nzl(k)=trapz(fy,squeeze(trapz(fx,squeeze(I_fft(k,:,:)),1)));
    x_beta(k)=-S_wavelength(num/2)*trapz(fy,squeeze(trapz(fx,FX.*squeeze(I_fft(k,:,:)),2)))/x_beta_nzl(k);
end
%计算Ax
for k=1:num
    x_d_Efft(k,:,:)=-i*2*pi*fftshift(fft2(squeeze(E_S_am(k,:,:)).*X)).*conj(squeeze(E_fft(k,:,:)));
    x_d_Efft(k,:,:)=-x_d_Efft(k,:,:)+conj(x_d_Efft(k,:,:));
    x_A(k)=i*S_wavelength(num/2)/(2*pi)*trapz(fy,squeeze(trapz(fx,FX.*squeeze(y_d_Efft(k,:,:)),1)))/x_beta_nzl(k);
end
%计算空间方差
for k=1:num
    x_theta(k)=trapz(y,squeeze(trapz(x,(X-x_av(k)).^2.*squeeze(I(k,:,:)),1)))/I_nlz(k);
end
%计算空间频率方差
for k=1:num
    fx_theta(k)=trapz(fx,(trapz(fy,(S_wavelength(num/2)*FX-x_beta(k)).^2.*squeeze(I_fft(k,:,:)),2)))/S_wavelength(num/2)^2/y_beta_nzl(k);
end
%计算Zy
Zx=(x_A+2*x_beta.*x_av)./(2*S_wavelength(num)^2*fx_theta);
x_theta_0=x_theta-Zx.^2*S_wavelength(num)^2.*fx_theta;
x_M_Squ=4*pi*x_theta_0.^0.5.*fx_theta.^0.5;
plot(zz(2:num+1)*1000,x_M_Squ)
hold on
toc