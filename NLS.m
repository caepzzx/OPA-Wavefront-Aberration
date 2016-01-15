% This code solves the NLS equation with the split-step method
% idu/dz-sgn(beta2)/2d^2u/d(tau)^2+N^2*|u|^2*u=0

%--Specify input parameters
clear all;%
distance=2;
beta2=-1;
N=1;%soliton order
mshape=3;
chirp0=0;%input pulse chirp(default value)

%---set simulation parameters
nt=1024;Tmax=32; %FFT points and window size
step_num=round(20*distance*N^2);
deltaz=distance/step_num; %step size in z
dtau=(2*Tmax)/nt; %step size in tau

%---tau and omega arrays
tau=(-nt/2:nt/2-1)*dtau; %temporal grid
omega=(pi./Tmax).*[(0:nt/2-1) (-nt/2:-1)];%frequency grid

%---Input Field profile
if mshape==0
    uu=sech(tau).*exp(-0.5i*chirp0*tau.^2); % soliton
else
    uu=exp(-0.5*(1+1i*chirp0).*tau.^(2*mshape));
end

%---Plot input pulse shape and spectrum
temp=fftshift(fft(uu)).*(nt*dtau)/sqrt(2*pi); % spectrum
figure; subplot(2,1,1);
    plot(tau,abs(uu).^2,'--k');hold on;
    axis([-20 20 0 inf]);
    xlabel('Normalized Time');
    ylabel('Normalized Power');
    title('Input and Output Pulse Shape and Spectrum');
    subplot(2,1,2);
    plot(fftshift(omega)./(2.*pi),(abs(temp)).^2,'k');hold on;
    axis([-5 5 0 inf]);
    xlabel('Normalied Frequency');
    ylabel('Spectral Power');

%--store dispersive phase shifts to speedup code
dispersion=exp(1i*0.5*beta2*omega.^2*deltaz);
hhz=1i*N^2*deltaz; %nonlinear phase factor

%*************[Beginning of MAIN LOOP]**************
%scheme:1/2N->D->1/2N ;first step nonlinear
temp=uu.*exp(abs(uu).^2*hhz/2); %note hhz/2
for n=1:step_num 
    f_temp=fft(temp).*dispersion;
    uu=ifft(f_temp);
    temp=uu.*exp(abs(uu).^2*hhz);
end
uu=temp.*exp(-abs(uu).^2*hhz/2);
temp=fftshift(fft(uu)).*(nt*dtau)/sqrt(2*pi); %Final spectrum

%**************[End of MAIN Loop]**********************
%---Plot output pulse shape and spectrum

subplot(2,1,1)
plot(tau,abs(uu).^2,'--r')
subplot(2,1,2)
plot(fftshift(omega)./(2.*pi),abs(temp).^2,'r')