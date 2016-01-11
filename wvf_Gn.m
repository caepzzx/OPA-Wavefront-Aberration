function Exy_ph=wvf_Gn(x,y,sgx)
%产生随机相位屏
%sgx 相位变化空间频率参数，通常2-12厘米

nx=size(x,2);
ny=size(y,2);
sgy=sgx;
[X,Y]=meshgrid(x,y);
a=2*rand(nx,ny)-1;
b=exp(-((X/sgx).^2+(Y/sgy).^2));
Exy_ph=conv2(a,b,'same');

% EFxy_ph=exp(-(X.^2+Y.^2)/(d0/5/log(2)^0.5)^2); 
% % subplot(2,2,1)
% % plot(x,EFxy_ph(nx/2,:))
% % hold on
% ns_p=unifrnd(-pi,pi,nx,ny);
% EFxy_ph=EFxy_ph.*exp(i*ns_p);
% [Exy_ph,fx,fy]=xy_fft(EFxy_ph,x,y);
% Exy_ph=abs(Exy_ph);
% Exy_ph=2*pi*(Exy_ph)/max(max(Exy_ph))-pi;
% 
% % figure(1)
% % subplot(2,1,1)
% % plot(x,Exy_ph(nx/2,:))
% % hold on
% % [phx,phy]=gradient(Exy_ph,dx);
% % subplot(2,1,2)
% % plot(x,phx(nx/2,:));
% % hold on