function Exy_ph=wvf_Gn(x,y,sgx,A)
%产生随机相位屏
%sgx 相位变化空间频率参数，通常2-12厘米
%A 控制产生相位的幅度

nx=size(x,2);
ny=size(y,2);
sgy=sgx;
[X,Y]=meshgrid(x,y);
a=2*rand(nx,ny)-1;
b=exp(-((X/sgx).^2+(Y/sgy).^2));
Exy_ph=A*conv2(a,b,'same');

