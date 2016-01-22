%
n=48;
I=E_S_out.*conj(E_S_out);
I_av_y=sum(squeeze(I(num/2,nx/2,96:150)))/n;
I_av_x=sum(squeeze(I(num/2,:,nx/2)))/49;
y_mod=(sum((squeeze((squeeze(I(num/2,nx/2,105:152))-I_av_y)).^2))/n/I_av_y^2).^0.5;
x_mod=(sum((squeeze((squeeze(I(num/2,96:144,nx/2))-I_av_x)).^2))/49/I_av_x^2).^0.5;