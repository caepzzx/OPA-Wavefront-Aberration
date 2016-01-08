function dyz=derivs(z,y,P_w,S_w,I_w,K_con,dk)
%[zzx]写出dy/dz表达式

n=size(y);
dyz=zeros(n);
dyz(1,:,:)=-1i*K_con(1).*y(3,:,:).*conj(y(2,:,:))*exp(-1i*dk*z);
dyz(2,:,:)=-1i*K_con(2).*y(3,:,:).*conj(y(1,:,:))*exp(-1i*dk*z);
dyz(3,:,:)=-1i*K_con(3).*y(1,:,:).*y(2,:,:)*exp(1i*dk*z);
