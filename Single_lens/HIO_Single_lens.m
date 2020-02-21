clc;
clear all;
close all;

folderpath = 'D:\Research_implementations\Paper_12\Single lens\Output\HIO';
a=dir([folderpath '/*.png']);
num=size(a,1)+1;

lambda=650*10^-9;
N=2^10;
Dimesnion=20*10^-3;
k=(2*pi)/lambda;
R1=Dimesnion/4;
R2=Dimesnion/5;
d1=Dimesnion/N;
alpha1=0.002;
beta1=0.0001;
alpha2=alpha1/1.1;
beta2=beta1/2;
cx1=0;
cy1=0;
cx2=-R1/2;
cy2=0*10^-3;
x=((-N/2):((N/2)-1))*d1;
y=((-N/2):((N/2)-1))*d1;
du = 1/(N*d1);
umax = 1/(2*d1);
u = -umax:du:umax-du;
[U,V]=meshgrid(u,u);

[X,Y]=meshgrid(x,y);
L1= real(sqrt(R1^2-((X-cx1).^2 + (Y-cy1).^2)));
L2= real(sqrt(R2^2-((X-cx2).^2 + (Y-cy2).^2)));
phase1 = k*alpha1*L1/1.4;
phase2 = k*alpha2*L2;
phase= phase1;
T1 = k*beta1*L1/8;
T2 = k*beta2*L2;

T= exp(-(T1));

fx=T;

figure
imagesc(x,y,T);
axis image
title(['Transmission profile of Object']);
colorbar
colormap(jet)

du = 1/(N*d1);
umax = 1/(2*d1);
u = -umax:du:umax-du;
[U,V]=meshgrid(u,u);
freq_sq=U.^2+V.^2;

noise=randn(N).*sqrt(10^-10);



mod_Fu=(abs(fftshift(fft2(fx))));


autocorrel=real(ifftshift(ifft2(((abs(fft2(fx)).^2)))));
figure 
imagesc(x,y,autocorrel);
title(['Autocorrelation-function']);
xlabel('x');
ylabel('y');
colorbar
colormap(jet)
% saveas(gcf,"3_autocorrelation.png");


M=zeros(N);
A=(X-cx1).^2+ (Y-cy1).^2<=(R1)^2;
%A= (abs(X)<=R1)&(abs(Y)<=R1 );

M(A)=1;
fx_estimate=rand(N).*M+imcomplement(M);

figure 
imagesc(x,y,fx_estimate);
title(['Initial estimate of f(x)']);
xlabel('x');
ylabel('y');
colorbar
colormap(jet)
%saveas(gcf,"4_fx_estimate.png");

it=30;
beta=0.1;

for i=1:it
F_estimate= fftshift(fft2(fx_estimate));
Ef(i)=(sum(sum(((abs(F_estimate)- mod_Fu).^2)*(du)*(du))))/(sum(sum((mod_Fu.^2)*(du)*(du))));
F_estimate=mod_Fu.*exp(1j.*angle(F_estimate));
fx_estimate_dash= real(ifft2(ifftshift(F_estimate)));ff=imcomplement(M).*fx_estimate_dash;
fff=ff~=1;
Es(i)=(sum(sum(((fx_estimate_dash(fff)).^2)*(d1)*(d1))))/(sum(sum((fx_estimate_dash.^2)*(d1)*(d1))));
dd=fx_estimate_dash>=0;
fx_estimate=(fx_estimate_dash.*dd + (fx_estimate-(beta.*fx_estimate_dash)).*imcomplement(dd) ).*M+imcomplement(M);
end

figure 
imagesc(x,y,(fx_estimate));
title(['Estimate of f(x)']);
xlabel('x');
ylabel('y');
colorbar
colormap(jet)
%saveas(gcf,"5_fx_estimate_final.png");


 i=1:it;
 
figure 
plot(i,Ef);
title(['Error E_F']);
xlabel('Iterations');
%saveas(gcf,"6_Ef.png");

figure 
plot(i,Es);
title(['Error E_s']);
xlabel('Iterations');
%saveas(gcf,"7_Es.png");
