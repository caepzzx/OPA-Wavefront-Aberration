function [Exy,fx,fy]=xy_fft(Exy,x,y)
% 本程序对输入二维数组进行傅立叶变换；
%-----------------------------
% Ey—为输入空域电场强度；
% y—为输入y轴坐标变量；
% Efy为对y方向进行空间傅立叶变换后的输出值；
% fy—为空间频率向量；
%-------------------------------
[nx,ny]=size(Exy);   
Y=y(ny)-y(1);      
dy=Y/ny;            %空间分辨率
dfy=1/Y;            %空间频率分辨率
fy=(-ny/2:(ny/2-1))*dfy;           %空间频率向量
X=x(nx)-x(1);      
dx=X/nx;            %空间分辨率
dfx=1/X;            %空间频率分辨率
fx=(-nx/2:(nx/2-1))*dfx;           %空间频率向量
Exy=fft2(Exy)*dx*dy; 
Exy=fftshift(Exy);  %把零频率移至中间位置