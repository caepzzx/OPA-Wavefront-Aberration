function Exy=xy_ifft(Exy,x,y)
% 本程序对输入二维数组的列向量进行傅立叶变换；
%-----------------------------
% Ey—为输入空域电场强度；
% y—为输入y轴坐标变量；
% Efy为对y方向进行空间傅立叶变换后的输出值；
% fy—为空间频率向量；
%-------------------------------
[nx,ny]=size(Exy);   
Y=y(ny)-y(1);      
dy=Y/(ny-1);            %空间分辨率
X=x(nx)-x(1);      
dx=X/(nx-1);            %空间分辨率
Exy=ifftshift(Exy);  %把零频率移至中间位置
Exy=ifft2(Exy)/dx/dy; 
Exy=ifftshift(Exy); 