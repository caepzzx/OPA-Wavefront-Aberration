function dyz=derivs(z,y,P_w,S_w,I_w,K_con,dk)
n=size(y);
dyz=zeros(n);
dyz(1,:,:)=-i*K_con.*y(3,:,:).*conj(y(2,:,:))*exp(-i*dk*z);
dyz(2,:,:)=-i*I_w/S_w*K_con.*y(3,:,:).*conj(y(1,:,:))*exp(-i*dk*z);
dyz(3,:,:)=-i*P_w/S_w*K_con.*y(1,:,:).*y(2,:,:)*exp(i*dk*z);
