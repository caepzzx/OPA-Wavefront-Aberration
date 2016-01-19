lambda=500e-9;     %设定波长，以Lambda表示波长
d=3e-6;            %设定光栅常数，以d表示光栅常数
a=1e-6;            %设定缝宽，以d表示缝宽
f=0.64;            %设定焦距，以f表示焦距
warning off MATLAB:divideByZero
N=4;              %设定缝数，以N表示缝数
ym=2*lambda*f/a;   %设定y方向的范围
xs=ym;            %设定x方向的范围
n=1001;
ys=linspace(-ym,ym,n);
%用线性采样法产生一个一维数组ys，n是此次采样的总点数 
%采样的范围从-ymax到ymax，采样的数组命名为ys 
%此数组装的是屏幕上的采样点的纵坐标 
for i=1:n           %对屏幕上的全部点进行循环计算，则要进行n次计算
  sinphi=ys(i)/f;  %以下几行进行光强的计算
  alpha=pi*a*sinphi/lambda;
  beta=pi*d*sinphi/lambda;
  B(i)=(sin(alpha)./alpha).^2.*(sin(N*beta)./sin(beta)).^2;
  B1=B/max(B);
end
NC=255;           %确定使用的灰度等级为255级
Br=(B/max(B))*NC;  %定标：使最大光强对应于最大灰度级（白色）
subplot(2,1,1);       %选中第一个子坐标轴
cla;                %清除轴上图形
image(xs,ys,Br);      %用image绘图函数创建图像
colormap(gray(NC));  %用灰度级颜色图设置色图和明暗
title('光栅衍射条纹');        %取名为光栅衍射条纹
subplot(2,1,2);       %选中第二个子坐标轴
cla;
plot(ys,B1);         %用plot函数绘制曲线
title('光栅衍射曲线');       %取名为光栅衍射曲线
