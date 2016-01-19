% function n=gain_smallsignal(lamdas0,lamdap0,material)

c=3e8;
pm=1e-12;
GW=1e9;
cm=1e-2;
v=1;
epi0=8.854e-12;
lamdap=.8e-6;
lamdas=1.2e-6;
material='BBO';
% lamdas=lamdas0;
% lamdap=lamdap0;
lamdai=1/(1/lamdap-1/lamdas);
no_s=sellmeier(lamdas,material,'o');
no_i=sellmeier(lamdai,material,'o');
no_p=sellmeier(lamdap,material,'o');
ne_p=sellmeier(lamdap,material,'e');

wp=2*pi*c/lamdap;
ws=2*pi*c/lamdas;
wi=2*pi*c/lamdai;

vgs=c/no_s;
vgi=c/no_i;

ks=no_s*ws/c;
ki=no_i*wi/c;

alpha=asind(sqrt((1-vgs^2/vgi^2)/(1+2*ks*vgs/(ki*vgi)+ks^2/ki^2)));
beta=acosd(vgs/vgi);
kp=sqrt(ks^2+ki^2-2*ks*ki*cosd(180-beta));
ne_pt=(kp*c/wp);
theta=asind((ne_p/ne_pt)*sqrt((no_p^2-ne_pt^2)/(no_p^2-ne_p^2)));
fprintf('泵浦光与光轴的夹角为:%f°\n',theta);
fprintf('泵浦光与信号光的夹角为:%f°\n',alpha);
fprintf('信号光与闲置光的夹角为:%f°\n',beta);

ni=no_i;
ns=no_s;
np=ne_pt;
kp=np*wp/c;
deff=2*pm/v;
Ip=linspace(10,100,100)*GW/cm^2;
gamma=sqrt(8*pi^2*deff^2*Ip/(ni*ns*np*lamdai*lamdas*epi0*c));
L=[1,2,3,4,5]*1e-3;
col='rgbcymk';
nn=length(gamma);
G=zeros(5,nn);
figure(1)
for i=1:5
G(i,:)=exp(2*gamma*L(i))/4;
semilogy(Ip/(GW/cm^2),G(i,:),'color',col(i),'LineWidth',2)
hold on
end
ylim([1,1e7])
legend('1mm','2mm','3mm','4mm','5mm',4)
xlabel('Pump Intensity(GW/cm^2)')
ylabel('Parametric Gain')
title(['lambda_p=',sprintf('%2.2f um, ',lamdap*1e6),'lambda_s=',sprintf('%2.2f um, ',lamdas*1e6),'material:',sprintf('%c',material)])

figure(2)
lamda=linspace(2e-7,4e-6,200);
no=sellmeier(lamda,material,'o');
ne=sellmeier(lamda,material,'e');
plot(lamda*1e6,abs(no),'r',lamda*1e6,abs(ne),'g','LineWidth',2)
legend('no','ne')
grid on
xlabel('wavelength/um')
ylabel('n')
title('BBO折射率曲线')

% figure(3)
% plot([0,0,1.1*kp*1e-6,1.1*kp*1e-6,0],[-5,5,5,-5,-5],'m','LineWidth',4)
% hold on
% arrow([0,0],[5*cosd(theta),5*sind(theta)],'tipangle',10,'width',1,'Length',10, 'EdgeColor','k','FaceColor','k')
% hold on
% arrow([0,0],[kp*1e-6,0],'tipangle',10,'width',1,'Length',10, 'EdgeColor','r','FaceColor','r')
% hold on
% arrow([0,0],[ks*1e-6*cosd(alpha),ks*1e-6*sind(alpha)],'tipangle',10,'width',1,'Length',10, 'EdgeColor','g','FaceColor','g')
% hold on
% arrow([kp*1e-6,0],[kp*(1e-6)+(ki*(1e-6)*cosd(180-beta+alpha)),ki*(1e-6)*sind(180-beta+alpha)],'tipangle',10,'width',1,'Length',10, 'EdgeColor','b','FaceColor','b')
% legend('晶体界面','光轴','Kp','Ks','Ki',5)
% title('相位匹配角')

%% 1共线的情况
% os+oi->ep
N=2000;%信号光的点数
M=2000;%角度的点数
lamdas=linspace(.3e-6,2.5e-6,N);
no_s=sellmeier(lamdas,material,'o');
ws=2*pi*c./lamdas;
lamdai=1./(1./lamdap-1./lamdas);
no_i=sellmeier(lamdai,material,'o');
wi=2*pi*c./lamdai;
wp=2*pi*c/lamdap;
no_p=sellmeier(lamdap,material,'o');
ne_p=sellmeier(lamdap,material,'e');
ne_pt=(no_s.*ws+no_i.*wi)./wp;
theta2=asind((ne_p./ne_pt).*sqrt((no_p.^2-ne_pt.^2)./(no_p.^2-ne_p.^2)));
vgs=c./no_s;
vgi=c./no_i;
vgp=c./ne_pt;
dip=1./vgi-1./vgp;
dsp=1./vgs-1./vgp;

figure(4)
plot(theta2,lamdas*1e6,'g','LineWidth',2)
hold on


%os+ei->ep
the=linspace(0,90,M);
no_s=sellmeier(lamdas,material,'o');
ws=2*pi*c./lamdas;

lamdai=1./(1./lamdap-1./lamdas);
no_i=sellmeier(lamdai,material,'o');
ne_i=sellmeier(lamdai,material,'e');
wi=2*pi*c./lamdai;

no_p=sellmeier(lamdap,material,'o');
ne_p=sellmeier(lamdap,material,'e');
wp=2*pi*c/lamdap;
ne_it0=zeros(1,N);
ne_pt0=zeros(1,N);
for i=1:N
ne_it=1./sqrt(sind(the).^2./ne_i(i)^2+cosd(the).^2./no_i(i).^2);
ne_pt=1./sqrt(sind(the).^2./ne_p^2+cosd(the).^2./no_p.^2);
[q,p]=min(abs(no_s(i)*ws(i)+ne_it*wi(i)-ne_pt*wp));
ne_it0(i)=ne_it(p);
ne_pt0(i)=ne_pt(p);
theta3(i)=the(p);
end

vgs=c./no_s;
vgi=c./ne_it0;
vgp=c./ne_pt0;
dip2=1./vgi-1./vgp;
dsp2=1./vgs-1./vgp;

plot(theta3,lamdas*1e6,'r','LineWidth',2)
hold on

% es+oi->ep
the=linspace(0,90,M);
no_s=sellmeier(lamdas,material,'o');
ne_s=sellmeier(lamdas,material,'e');
ws=2*pi*c./lamdas;

lamdai=1./(1./lamdap-1./lamdas);
no_i=sellmeier(lamdai,material,'o');
wi=2*pi*c./lamdai;

no_p=sellmeier(lamdap,material,'o');
ne_p=sellmeier(lamdap,material,'e');
wp=2*pi*c/lamdap;
ne_st00=zeros(1,N);
ne_pt00=zeros(1,N);

for i=1:N
ne_st=1./sqrt(sind(the).^2./ne_s(i)^2+cosd(the).^2./no_s(i).^2);
ne_pt=1./sqrt(sind(the).^2./ne_p^2+cosd(the).^2./no_p.^2);
[q,p]=min(abs(ne_st*ws(i)+no_i(i)*wi(i)-ne_pt*wp));
ne_st00(i)=ne_st(p);
ne_pt00(i)=ne_pt(p);
theta3(i)=the(p);
end
vgs=c./ne_st00;
vgi=c./no_i;
vgp=c./ne_pt00;
dip3=1./vgi-1./vgp;
dsp3=1./vgs-1./vgp;
plot(theta3,lamdas*1e6,'b','LineWidth',2)
xlim([18,40])
ylim([1,2.5])
legend('os+oi->ep','os+ei->ep','es+oi->ep')

figure(5)
plot(lamdas*1e6,dip*1e12,lamdas*1e6,dsp*1e12,'r');
hold on
plot(lamdas*1e6,dip2*1e12,lamdas*1e6,dsp2*1e12,'g');
hold on
plot(lamdas*1e6,dip3*1e12,lamdas*1e6,dsp3*1e12,'b');
xlim([1,1.6])
ylim([-150,50])
% end

