function [bmx,bmy]=Beam_Quality(x,y,fx,fy,wlgth,Exy)
% Beam quality calculation
% x,y                  横向坐标
% Exy                  输入场
% betax,betay          x,y向光束倾斜角（频域内）
% cx,cy                质心坐标
% thetax,thetay        x,y向二阶矩(开平方）
% thetafx,thetafy      fx,fy向二阶矩(开平方）
% bmx,bmy              光束质量因子
% 编写于2008/12/10
[nt,nx,ny]=size(Exy);
[X,Y]=meshgrid(x,y);
[FX,FY]=meshgrid(fx,fy);
bmx=zeros(1,nt);
bmy=zeros(1,nt);
cx=zeros(1,nt);
cy=zeros(1,nt);
betax=zeros(1,nt);
betay=zeros(1,nt);
thetax=zeros(1,nt);
thetay=zeros(1,nt);
thetafx=zeros(1,nt);
thetafy=zeros(1,nt);
Ifft=zeros(nt,nx,ny);
Efft=zeros(nt,nx,ny);
y_d_Efft=zeros(nt,nx,ny);
x_d_Efft=zeros(nt,nx,ny);
y_A=zeros(1,nt);
x_A=zeros(1,nt);
If=zeros(1,nt);

I=Exy.*conj(Exy);
It=trapz(x,trapz(y,I,3),2);
for k=1:nt
    cy(k)=trapz(x,trapz(y,Y.*squeeze(I(k,:,:)),2))/It(k);                   %质心坐标
    Efft(k,:,:)=fftshift(fft2(squeeze(Exy(k,:,:))));
    Ifft(k,:,:)=Efft(k,:,:).*conj(Efft(k,:,:));
    If(k)=trapz(fx,trapz(fy,squeeze(Ifft(k,:,:)),2));
    betay(k)=-wlgth*trapz(fx,trapz(fy,FY.*squeeze(Ifft(k,:,:)),2))/If(k);
    
    y_d_Efft(k,:,:)=i*2*pi*conj(fftshift(fft2(squeeze(Exy(k,:,:)).*Y))).*squeeze(Efft(k,:,:));%计算Ay
    y_d_Efft(k,:,:)=y_d_Efft(k,:,:)-conj(y_d_Efft(k,:,:));
    y_A(k)=i*wlgth/(2*pi)*trapz(fx,squeeze(trapz(fy,FY.*squeeze(y_d_Efft(k,:,:)),2)))/If(k);
    
    thetay(k)=trapz(x,trapz(y,(Y-cy(k)).^2.*squeeze(I(k,:,:)),2))/It(k);
    thetafy(k)=trapz(fx,trapz(fy,(wlgth*FY-betay(k)).^2.*squeeze(Ifft(k,:,:)),2))/wlgth^2/If(k);
    
    cx(k)=trapz(x,trapz(y,X.*squeeze(I(k,:,:)),2))/It(k);                   %质心坐标
    Efft(k,:,:)=fftshift(fft2(squeeze(Exy(k,:,:))));
    Ifft(k,:,:)=Efft(k,:,:).*conj(Efft(k,:,:));
    betax(k)=-wlgth*trapz(fx,trapz(fy,FX.*squeeze(Ifft(k,:,:)),2))/If(k);
    
    x_d_Efft(k,:,:)=i*2*pi*conj(fftshift(fft2(squeeze(Exy(k,:,:)).*X))).*squeeze(Efft(k,:,:));%计算Ay
    x_d_Efft(k,:,:)=y_d_Efft(k,:,:)-conj(y_d_Efft(k,:,:));
    x_A(k)=i*wlgth/(2*pi)*trapz(fx,squeeze(trapz(fy,FX.*squeeze(y_d_Efft(k,:,:)),2)))/If(k);
    
    thetax(k)=trapz(x,trapz(y,(X-cx(k)).^2.*squeeze(I(k,:,:)),2))/It(k);
    thetafx(k)=trapz(fx,trapz(fy,(wlgth*FX-betax(k)).^2.*squeeze(Ifft(k,:,:)),2))/wlgth^2/If(k);
end
Zy=(y_A+2*betay.*cy)./(2*wlgth^2*thetafy);
thetay0=thetay-Zy.^2*wlgth^2.*thetafy;
bmy=4*pi*thetay0.^0.5.*thetafy.^0.5;
%计算x向光束质量

Zx=(x_A+2*betax.*cx)./(2*wlgth^2*thetafx);
thetax0=thetax-Zx.^2*wlgth^2.*thetafx;
bmx=4*pi*thetax0.^0.5.*thetafx.^0.5;

