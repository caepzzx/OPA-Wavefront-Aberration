function Exy_ph=wvf_Gn(x,y,d0,dx)
nx=size(x,2);
ny=size(y,2);
[X,Y]=meshgrid(x,y);
EFxy_ph=exp(-(X.^2+Y.^2)/(d0/5/log(2)^0.5)^2); 
% subplot(2,2,1)
% plot(x,EFxy_ph(nx/2,:))
% hold on
ns_p=unifrnd(-pi,pi,nx,ny);
EFxy_ph=EFxy_ph.*exp(i*ns_p);
[Exy_ph,fx,fy]=xy_fft(EFxy_ph,x,y);
Exy_ph=abs(Exy_ph);
Exy_ph=2*pi*(Exy_ph)/max(max(Exy_ph))-pi;

% figure(1)
% subplot(2,1,1)
% plot(x,Exy_ph(nx/2,:))
% hold on
% [phx,phy]=gradient(Exy_ph,dx);
% subplot(2,1,2)
% plot(x,phx(nx/2,:));
% hold on