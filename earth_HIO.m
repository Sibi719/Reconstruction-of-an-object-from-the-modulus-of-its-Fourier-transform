clc;
clear all;
close all;
N=1024;

Dimesnion=20*10^-3;
R1=Dimesnion/5.7;


fx=imread("4.jpg");
fx=rgb2gray(fx);
Nr=size(fx,1);
Nc=size(fx,2);
Dr=(N-Nr)/2;
Dc=(N-Nc)/2;
fx = padarray(fx,[Dr,Dc],0,'both');


d1=Dimesnion/N;
x=((-N/2):((N/2)-1))*d1;
y=((-N/2):((N/2)-1))*d1;
du = 1/(N*d1);
umax = 1/(2*d1);
u = -umax:du:umax-du;
[U,V]=meshgrid(u,u);
[X,Y]=meshgrid(x,y);

figure
imagesc(fx);
colormap(gray);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) ;
ax_height = outerpos(4) ;
ax.Position = [0 0 ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(gcf,"1_Input object.png");


folderpath = 'D:\Research_implementations\Paper_12\Earth_1\output\HIO';
a=dir([folderpath '/*.png']);
num=size(a,1)+1;

noise=randn(N).*sqrt(10^-10);

mod_Fu=(abs(fftshift(fft2(fx))));

M=zeros(N);
A=X.^2+ Y.^2<=(R1)^2;
%A= (abs(X)<=R1/2)&(abs(Y)<=R1/2 );
M(A)=1;

fx_estimate=256.*rand(N).*M ;
%save iniest.mat fx_estimate;

figure 
imagesc(x,y,fx_estimate);
colormap(gray);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) ;
ax_height = outerpos(4) ;
ax.Position = [0 0 ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(gcf,"3_fx_estimate.png");

it=2000;
beta=0.7;
for i=1:it
F_estimate= fftshift(fft2(fx_estimate));
Ef(i)=(sum(sum(((abs(F_estimate)- mod_Fu).^2)*(du)*(du))))/(sum(sum((mod_Fu.^2)*(du)*(du))));
F_estimate=mod_Fu.*exp(1j.*angle(F_estimate));
fx_estimate_dash= real(ifft2(ifftshift(F_estimate)));
Es(i)=(sum(sum(((fx_estimate_dash(fx_estimate_dash<0)).^2)*(d1)*(d1))))/(sum(sum((fx_estimate_dash.^2)*(d1)*(d1))));
dd=fx_estimate_dash>=0;
fx_estimate=(fx_estimate_dash.*dd)+((fx_estimate-(beta.*fx_estimate_dash)).*imcomplement(dd));
end

figure 
imagesc(fx_estimate);
colormap(gray)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) ;
ax_height = outerpos(4) ;
ax.Position = [0 0 ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(gcf,"output\HIO\"+num+".png");

i=1:it;
 
figure 
plot(i,Ef);
title(['Error E_F']);
xlabel('Iterations');
saveas(gcf,"4_Ef.png");

figure 
plot(i,Es);
title(['Error E_s']);
xlabel('Iterations');
saveas(gcf,"5_Es.png");

autocorrel=(ifft2((abs(fft2(fx)).^2)));

% figure 
% imagesc(x,y,autocorrel);
% title(['Autocorrelation-function']);
% xlabel('x');
% ylabel('y');
% colorbar
% colormap(gray)
% saveas(gcf,"2_autocorrelation.png");
% 