function [fx,fy]=spati_vector(nx,ny,x,y)
%产生空间频率向量 
Y=y(ny)-y(1);      
dfy=1/Y;            %空间频率分辨率
fy=(-ny/2:(ny/2-1))*dfy;           %空间频率向量
X=x(nx)-x(1);      
dfx=1/X;            %空间频率分辨率
fx=(-nx/2:(nx/2-1))*dfx;           %空间频率向量
end