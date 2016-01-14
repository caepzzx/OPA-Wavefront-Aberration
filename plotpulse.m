figure;
I=(abs(E_S_out(:,nx/2,ny/2))).^2*(1/2*c*ele_c).*S_R_index(:);
plot(t*1e9,I/max(I),'b','LineWidth',1);
hold on
I=(abs(E_I_out(:,nx/2,ny/2))).^2*(1/2*c*ele_c).*I_R_index(:);
plot(t*1e9,0.9*I/max(I),'g','LineWidth',1);
I=(abs(E_P_out(:,nx/2,ny/2))).^2*(1/2*c*ele_c).*P_R_index(:);
plot(t*1e9,I/max(I),'r','LineWidth',1);