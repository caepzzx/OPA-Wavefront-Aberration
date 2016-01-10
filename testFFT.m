wx=0.1;
wy=0.5;
L=2;
M=200;
% x=linspace(-L/2,L/2,M);  %x£­×ø±ê
dx=L/M;
x=-L/2:dx:L/2-dx;
y=x;
[X,Y]=meshgrid(x,y);
g=rect(X/(2*wx)).*rect(Y/(2*wy)); %Signal
figure(1)
imagesc(x,y,g);
colormap('gray');
axis square;
axis xy
xlabel('x (m)');
ylabel('y (m)');
[G,fx,fy]=xy_fft(g,x,y);
figure(2)
surf(fx,fy,abs(G))
camlight left; lighting phong
colormap('gray')
shading interp
ylabel('fy £¨cyc/m£©'); xlabel('fx (cyc/m)');
figure(3)
imagesc(x,y,xy_ifft(G,x,y));
colormap('gray');
axis square;
axis xy
xlabel('x (m)');
ylabel('y (m)');