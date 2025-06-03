%clc; clear; 
close all

pA=2.4882;EI=0.236;g=9.81;
l=0.6;
mw=13.1;
mk=0.2;
tmax=15;
% r=5000;
r=5000;
tmax = 15;
%r = 1333333;
dt=tmax/(r-1);


% cac buoc sai phan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                 % n=3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=9;
ds=l/(n-1); 
S1=(EI*(dt^2))/(mw*(ds^3));
S2=(-EI*(dt^2))/((ds^4)*pA);
S3=(EI*(dt^2))/(mk*ds^3);
F1=zeros(1,r);
F1(1:500)=10;
F1(501:1000)=-10;
w=zeros(n,r);
% mo phong diem x1
A = [-0.5 30 -0.00005 0.03];
% B = [100, 100, 320, 200];

for j =1:r-2
wsss0=(w(3,j+1)-2*w(2,j+1)+w(1,j+1))/(2*ds^3);

w(1,j+2)=(F1(j+2)-EI*wsss0)*(dt^2/mw)+2*w(1,j+1)-w(1,j);

% w(2,j+2)=(w(1,j+2)+w(3,j+2))/2;
for i=3:n-2
S2=(-EI*dt^2)/(ds^4*pA);
wssss=w(i+2,j+1)-4*w(i+1,j+1)+6*w(i,j+1)-4*w(i-1,j+1)+w(i-2,j+1);

w(i,j+2)=S2*wssss+2*w(i,j+1)-w(i,j);

end
w(2,j+2)=(w(1,j+2)+w(3,j+2))/2;
wsssl=(1*w(n,j+1)-3*w(n-1,j+1)+2*w(n-2,j+1));
w(n,j+2)=2*w(n,j+1)-w(n,j)+S3*wsssl;
w(n-1,j+2)=(w(n,j+2)+w(n-2,j+2))/2;
end
w = w;
% w = w/6;

x=0:dt:tmax;
%axes('position',[.8    .08    .5    .5]);
B = [100, 100, 600, 300];
figure('Position', B);
figure(1) %vi tri xe con
grid on 
hold on
plot(x,w(1,:), 'color',[0.07,0.62,1.00],'linewidth',2);
% legend(num2str(n),num2str(n1), num2str(n3), num2str(n4));
%title('Position of driving unit(m)')
xlabel('Thời gian (s)');
ylabel('Vị trí xe con (m)');
axis([0 10 0 2]);
xticks([0 2 4 6 8 10]);
yticks([0 0.5 1 1.5 2]);

figure('Position', B);
figure(2)%diem giua
hold on 
grid on
%plot(x,w(ceil(n/2),:)-w(1,:),'r', linewidth = 0.7);
% plot(x,w1(ceil(n1/2),:)-w1(1,:),'r');
% plot(x,w3(ceil(n3/2),:)-w3(1,:),'y');
% plot(x,w4(ceil(n4/2),:)-w4(1,:),'g');
% legend(num2str(n),num2str(n1), num2str(n3), num2str(n4));
%title('Deflection of the stacker crane');
% xlabel('Time(s)');
% C = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
% B = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
% plot(C,B);
% figure(3)%diem cuoi
% hold on 
% grid on
plot(x,w(n,:)-w(1,:), 'color',[0.07,0.62,1.00],'linewidth',2);
%plot(x,w(5,:)-w(1,:), linewidth = 1.2);
% %plot(x,w(n,:)-w(1,:),'b-.');
% plot(x,w1(n1,:)-w1(1,:),'r');
% plot(x,w3(n3,:)-w3(1,:),'m');
% plot(x,w4(n4,:)-w4(1,:),'g');
% labels = [n n1 n3 n4];
% legend(labels)
% legend(num2str(n),num2str(n1), num2str(n3), num2str(n4));
%title('Deflection of the end point of the stacker crane')
xlabel('Thời gian (s)');
ylabel('ω(L,t) (m)');
axis([0 10 -0.09 0.09]);
xticks([0 2 4 6 8 10]);
yticks([-0.09 -0.06 -0.03 0 0.03 0.06 0.09]);
legend('Đỉnh thang');


figure(4)
Fs=1/dt;
time=0:dt:(r-2)*dt-dt;
l1=r;
% subplot(2,2,1)
%fft_w=fft(out.W2,l1)*(2/l1);
fft_w=fft(w(n,:)-w(1,:),l1)*(2/l1);
abs_w=abs(fft_w);
freq=0:(1/time(end)):Fs/2-(1/time(end));
plot(freq,abs_w(1:length(freq)),'LineWidth',0.8); 
grid on
xlabel('Frequency(Hz)');
ylabel('Amplitude(m)');
title('Frrequency spectrum of vibration');
% axis(A);
% 
% 
% % 
figure('Position', [100 100 600 400]);
figure(3);
xc = w-w(1,:);
xc1 = xc';
u0 =linspace(0,l,n); t_tr=linspace(0,tmax,r);
surfc(u0,t_tr,xc1); view(65,10);
title({'Displacement of the beam'});
ylabel('t(s)'); xlabel('Y(m)', 'Fontsize', 12); zlabel('ω(Y,t)(m)','Fontsize', 12);
% axis([0 0.6 0 10 -0.02 0.02]);

figure('Position', [100 100 600 400]);
[X, Y] = meshgrid(0:dt:tmax,0:ds:l);
meshc(Y, X, xc);
ylabel('t(s)');
xlabel('Y(m)');
zlabel('ω(Y,t)(m)');
view(65,10);
axis([0 0.6 0 10 -0.08 0.08]);
xticks([0 0.2 0.4 0.6]);
yticks([0 2 4 6 8 10]);
