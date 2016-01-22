function [E_S_am,E_I_am,E_P_am,E_S_ph,E_I_ph,E_P_ph,sum_Ep,sum_Es]=Recod_data(E_S_out,E_I_out,E_P_out,num,x,y,t,t0)
% 记录相位，能量随传播距离的变化
% k 距离z的序数
% num 时间取样个数
% E_S_out，E_I_out，E_P_out 序数k对应距离处三波包络的电场
% E_S_am,E_I_am,E_P_am  序数k对应距离处三波包络的电场在中间时间序数的复振幅
% E_S_ph，E_I_ph，E_P_ph 序数k对应距离处三波的相位
% sum_Ep，sum_Es 序数k对应距离处p光和s光的能量
    const_LBO;
    
    E_S_am=E_S_out(num/2,:,:);
    E_I_am=E_I_out(num/2,:,:);
    E_P_am=E_P_out(num/2,:,:);
    E_S_ph=atan2(imag(E_S_out(num/2,:,:)),real(E_S_out(num/2,:,:)));
    E_I_ph=atan2(imag(E_I_out(num/2,:,:)),real(E_I_out(num/2,:,:)));
    E_P_ph=atan2(imag(E_P_out(num/2,:,:)),real(E_P_out(num/2,:,:)));
    %记录信号光能量随参量作用过程的变化
    Is=(1/2*c*S_R_index(num/2)*ele_c).*E_S_out.*conj(E_S_out);
    Ip=(1/2*c*P_R_index*ele_c).*E_P_out.*conj(E_P_out);
    sum_Ep=trapz(y,squeeze(trapz(x,squeeze(trapz(t,Ip,1)),1)))*1000;
    sum_Es=trapz(y,squeeze(trapz(x,squeeze(trapz(t,Is,1)),1)))*1000;
    