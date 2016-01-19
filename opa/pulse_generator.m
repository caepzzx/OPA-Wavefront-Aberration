function [At,intense] = pulse_generator(energy,diameter,type,tao_p,n,chirp0,t,m)
%PULSE_GENERATOR 此处显示有关此函数的摘要
%   此处显示详细说明

c=3e8;
epislon=8.85e-12;

switch type
    case 'gauss'
        A0=1;
        dt=t(2)-t(1);
        t0=tao_p./(2*sqrt(log(2)));
        At=A0*exp(-(1+1i*chirp0)*t.^2/(2*t0^2));
        norm=sum(0.5*epislon*c*n*abs(At).^2*dt);
        area=pi*diameter^2/4;
        At=At*sqrt(energy/(area*norm));
        intense=0.5*epislon*c*n*abs(At).^2;
        
    case 'supergauss'
        A0=1;
        t0=tao_p./(2*(log(2)).^(1./(2*m)));
        At=A0*exp(-(1+1i*chirp0)*((t./t0).^(2*m))/2);
        dt=t(2)-t(1);
        norm=sum(0.5*epislon*c*n*abs(At).^2*dt);
        area=pi*diameter^2/4;
        At=At*sqrt(energy/(area*norm));      
        intense=0.5*epislon*c*n*abs(At).^2;

        
    case 'sech'
        A0=1;
        t0=tao_p./(2*log(1+sqrt(2)));
        At=A0*sech(t./t0).*exp(-1i*chirp0*t.^2/(2*t0^2));
        dt=t(2)-t(1);
        norm=sum(0.5*epislon*c*n*abs(At).^2*dt);
        area=pi*diameter^2/4;
        At=At*sqrt(energy/(area*norm));       
        intense=0.5*epislon*c*n*abs(At).^2;
        
end

