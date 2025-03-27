clc;
clear;
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pA=0.297;EI=0.754;g=9.81;
L=0.7;mw=13.1;mh=0.86;mk=0.04;m=14.15;
c=0;
K_0=1;%gain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chọn n,r, thời gian mô phỏng
tmax=15;
n=9;r=15000;
%%%%%%%%%%%%%%%%%%%%%%%%
dt=tmax/(r-1);
ds=L/(n-1); 
e=zeros(1,r);
F1=zeros(1,r);
F1_is=zeros(1,r);
h=zeros(1,r);h(1:3)=ds;
h_is=zeros(1,r);h_is(1:r)=ds;
w3=zeros(1,r);
w3_is=zeros(1,r);
F2=zeros(1,r);
l=2;%% vị trí ban đầu của xe nâng
l_is=2;
w=zeros(n,r);
w_is=zeros(n,r);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chọn bộ điều khiển + SP 
% 0 - PD 
% 1 - PD IS tần thay đổi
% 2 - PD IS thường
sp1=1;
sp2=0.5;
flag=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Chọn bộ thông số IS
%ETM4
% f=6.13415;
sp1_VTIS=zeros(1,r);
w_0=6.13415*2*3.14;
% w_0=4.66729*2*3.14;
f=w_0/(2*3.14);
k_is=exp(c*3.14/sqrt(1-c^2));
m_is=1;
i_is=(1+m_is)*k_is^2/(k_is^2+(1+m_is)*(k_is^(4/3)+k_is^(2/3))+m_is);
a1=i_is/(1+m_is);
a2=i_is/k_is^(2/3);
a3=i_is/k_is^(4/3);
a4=m_is*i_is/((1+m_is)*k_is^2);
t1=0;
t2=(2*3.14/3)/(2*3.14*f);
t3=(4*3.14/3)/(2*3.14*f);
t4=(2*3.14)/(2*3.14*f);
kp1 =8; kd1 =10 ;
kp2 =20; kd2 =30;
canc=sqrt(1-c^2);
K_t=K_0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       MÔ PHỎNG
% for k =3:r-1
% wsss0_is=(w_is(3,k)-2*w_is(2,k)+w_is(1,k))/(2*ds^3);
% w_is(1,k+1)=(F1_is(k+1)-EI*wsss0_is)*(dt^2/mw)+2*w_is(1,k)-w_is(1,k-1);
% for m=3:n-2
% S2=(-EI*dt^2)/(ds^4*pA);
% S21=(-EI*dt^2)/ds^3*mh;
% C=-c*dt^2*(w_is(m,k)-w_is(m,k-1))/dt;
% wssss_is=w_is(m+2,k)-4*w_is(m+1,k)+6*w_is(m,k)-4*w_is(m-1,k)+w_is(m-2,k);
% ddtx3_is=(w3_is(k)-2*w3_is(k-1)+w3_is(k-2))/(dt^2);
% dyx2_is=(w_is(l+1,k-1)-w_is(l-1,k-1))/2*ds;
% h_is(k+1)=ds;
% l_is=2;
% w_is(m,k+1)=S2*wssss_is+C+2*w_is(m,k)-w_is(m,k-1);
% w3_is(k+1)=0;
%  S3=(EI*(dt^2))/(mk*2*ds^3);
%  w_is(2,k+1)=w_is(1,k+1);
%  wsssl_is=(-2*w_is(n,k)+3*w_is(n-1,k)-w_is(n-2,k));
%  w_is(n,k+1)=2*w_is(n,k)-w_is(n,k-1)+S3*wsssl_is;
%  w_is(n-1,k+1)=(w_is(n,k+1)+w_is(n-2,k+1))/2;
% sp1_0is=sp1*(a1*us(dt*k-t1)+a2*us(dt*k-t2)+a3*us(dt*k-t3)+a4*us(dt*k-t4));
% F1_is(k+2)= kp1*(sp1_0is-w_is(1,k+1)) + kd1*(w_is(1,k)-w_is(1,k+1))/dt;
% end
% end
for j =3:r-1
wsss0=(w(3,j)-2*w(2,j)+w(1,j))/(2*ds^3);
w(1,j+1)=(F1(j+1)-EI*wsss0)*(dt^2/mw)+2*w(1,j)-w(1,j-1);
% w(1,j+1)=(F1(j+1)-EI*wsss0+2*sin(j*pi/4))*(dt^2/mw)+2*w(1,j)-w(1,j-1);

for i=3:n-2
S2=(-EI*dt^2)/(ds^4*pA);
S21=(-EI*dt^2)/ds^3*mh;
C=-c*dt^2*(w(i,j)-w(i,j-1))/dt;
wssss=w(i+2,j)-4*w(i+1,j)+6*w(i,j)-4*w(i-1,j)+w(i-2,j);
ddtx3=(w3(j)-2*w3(j-1)+w3(j-2))/(dt^2);
dyx2=(w(l+1,j-1)-w(l-1,j-1))/2*ds;
h(j+1)=(F2(j+1)-mh*g-mh*ddtx3*dyx2)*dt^2/mh+2*h(j)-h(j-1);
if h(j+1)<ds
    l=2;
    h(j+1)=ds;
end
if h(j+1)>L-ds
    l=n-1;
    h(j+1)=L-ds;
end
if h(j+1)> ds && h(j+1)<L-ds
    l= ceil(h(j+1)/ds);
end
if i == l && l > 2 
w(i,j+1) =S21*wssss+C+2*w(i,j)-w(i,j-1);
w3(j+1)=w(i,j+1)-w(1,j+1);
end
if i ~= l && l >2 
w(i,j+1)=S2*wssss+C+2*w(i,j)-w(i,j-1);
w3(j+1)=w(i,j+1)-w(1,j+1);
end
if l==2
w(i,j+1)=S2*wssss+C+2*w(i,j)-w(i,j-1);
w3(j+1)=0;
end
 S3=(EI*(dt^2))/(mk*2*ds^3);
 w(2,j+1)=w(1,j+1);
 wsssl=(-2*w(n,j)+3*w(n-1,j)-w(n-2,j));
 w(n,j+1)=2*w(n,j)-w(n,j-1)+S3*wsssl;
 w(n-1,j+1)=(w(n,j+1)+w(n-2,j+1))/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag ==0
   % PD
   % P - đáp ứng nhanh, gây ra quá điều chỉnh
   % I - Giảm sai lệch tĩnh, gây ra quá điều chỉnh
   % D - Giảm thời gian xác lập, Dễ mất ổn định
% thông số
% kp1 =2.5; kd1 =9 ;
% kp2 =10; kd2 =20;
kp1=15;kd1=20;
kp2=15;kd2=20;
F1(j+2)= kp1*(sp1-w(1,j+1)) + kd1*(w(1,j)-w(1,j+1))/dt;
F2(j+2)=mh*g+ kp2*(sp2-h(j+1)) + kd2*(h(j)-h(j+1))/dt;
end
%PD+vary-time IS
if flag ==1
kp1 =10; kd1 =20 ;
% kp1 =8; kd1 =10 ;
kp2 =20; kd2 =30;
canc=sqrt(1-c^2);
K_t=K_0;
sp1_00is=sp1*(a1*us(dt*j-t1)+a2*us(dt*j-t2)+a3*us(dt*j-t3)+a4*us(dt*j-t4));
    if h(j+1)<3*ds
    w_t=w_0;
end
if h(j+1)>=3*ds && h(j+1)<4*ds
    w_t=4.66729*2*3.14;
end
if h(j+1)>=4*ds && h(j+1)<5*ds
    w_t=2.80037*2*3.14;
end
if h(j+1)>=5*ds && h(j+1)<6*ds
    w_t=2.06694*2*3.14;
end
if h(j+1)>=6*ds && h(j+1)<7*ds
    w_t=1.86692*2*3.14;
end
if h(j+1)>=7*ds
    w_t=2.06694*2*3.14;
end
% A_VTIS=(w_t^2-w_0^2)*(a1*K_0*(1-exp(-c*w_0*(dt*j-t1))/canc*sin(w_0*canc*(dt*j-t1)+acos(c)))+a2*K_0*(1-exp(-c*w_0*(dt*j-t2))/canc*sin(w_0*canc*(dt*j-t2)+acos(c)))+a3*K_0*(1-exp(-c*w_0*(dt*j-t3))/canc*sin(w_0*canc*(dt*j-t3)+acos(c)))+a4*K_0*(1-exp(-c*w_0*(dt*j-t4))/canc*sin(w_0*canc*(dt*j-t4)+acos(c))));
% B_VTIS=2*(c*w_t-c*w_0)*((a1*K_0*exp(-c*w_0*(dt*j-t1))/canc*sin(w_0*canc*(dt*j-t1)))+(a2*K_0*exp(-c*w_0*(dt*j-t2))/canc*sin(w_0*canc*(dt*j-t2)))+(a3*K_0*exp(-c*w_0*(dt*j-t3))/canc*sin(w_0*canc*(dt*j-t3)))+(a4*K_0*exp(-c*w_0*(dt*j-t4))/canc*sin(w_0*canc*(dt*j-t4))));
% C_VTIS=K_0*w_0^2*(a1*us(dt*j-t1)+a2*us(dt*j-t2)+a3*us(dt*j-t3)+a4*us(dt*j-t4));
if (j+1)==r
    break
end
F1(j+2)=kp1*(sp1_00is-w(1,j+1)) + kd1*(w(1,j)-w(1,j+1))/dt;
wsss0=(w(3,j+1)-2*w(2,j+1)+w(1,j+1))/(2*ds^3);
w(1,j+2)=(F1(j+2)-EI*wsss0)*(dt^2/mw)+2*w(1,j+1)-w(1,j);
A_VTIS=(w_t^2-w_0^2)*w(1,j+2);
B_VTIS=2*c*(w_t-w_0)*(w(1,j+2)-w(1,j))/dt;
C_VTIS=K_0*w_0^2*sp1_00is;
sp1_VTIS=1/(K_t*w_t^2)*(A_VTIS+B_VTIS+C_VTIS);
F1(j+2)= kp1*(sp1_VTIS-w(1,j+1)) + kd1*(w(1,j)-w(1,j+1))/dt;
% F2(j+2)=mh*g+ kp2*(sp2-h(j+1)) + kd2*(h(j)-h(j+1))/dt; 
end
if flag ==2
%PD+IS
ki1 =8; kd1 =10 ;
% ki2 =1700; kd2 =1200;
ki2 =10; kd2 =20;
sp1IS=sp1*(a1*us(dt*j-t1)+a2*us(dt*j-t2)+a3*us(dt*j-t3)+a4*us(dt*j-t4));
F1(j+2)= ki1*(sp1IS-w(1,j+1)) + kd1*(w(1,j)-w(1,j+1))/dt;
% F2(j+2)=mh*g+ ki2*(sp2-h(j+1)) + kd2*(h(j)-h(j+1))/dt; 
end
end
end
x=0:dt:tmax;
figure(1)
subplot(2,2,1)
grid on 
hold on
plot(x,w(1,:),'b',LineWidth=2.5);
xlabel('Thời gian(s)')
ylim([0 2])
ylabel('Vị trí xe con(m)')
title('Vị trí xe con')
subplot(2,2,2)
hold on 
grid on
dolac=zeros(1,r);
for i=1:r
dolac(i)=w(n,i)-w(1,i);
end
plot(x,dolac,'b',LineWidth=1);
title('Độ lắc điểm cuối');
xlabel('Thời gian(s)');
ylabel('m');
axis([0 10 -0.01 0.01])
check = w3-w(1,:);
subplot(2,2,3)
plot(x,h,'r');
title('Vị trí xe nâng',LineWidth=1);
xlabel('Thời gian(s)');
ylabel('m');
title('Vị trí xe nâng');
grid on
hold on
subplot(2,2,4)
grid on
hold on
plot(x,w3,LineWidth=1)
title('Độ lắc xe nâng');
xlabel('Thời gian(s)');
ylabel('m');
figure(2)
subplot(1,2,1)
plot(x,F1(1:r),LineWidth=1);
title('F1')
grid on
hold on
subplot(1,2,2)
plot(x,F2(1:r),LineWidth=1);
title('F2');
grid on
hold on

% Phân tích phổ
% figure(3)
% for i=1:n
% Fs=1/dt;
% time=0:dt:(r-2)*dt-dt;
% l1=r;
% fft_w=fft(w(i,:)-w(1,:),l1)*(2/l1);
% abs_w=abs(fft_w);
% freq=0:(1/time(end)):Fs/2-(1/time(end));
% subplot(ceil(n/2), 2, i);
% plot(freq,abs_w(1:length(freq)),'LineWidth',0.8);
% hold on
% grid on
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 9);
% xlabel('Frequency(Hz)');
% ylabel('Amplitude(m)');
% title('Pho dao dong');
% title(['FFT tại vị trí x = ', num2str(ds*(i-1)), ' m']);
% end
% 
