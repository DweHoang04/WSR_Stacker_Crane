clc;clear;close all
pA=0.297;EI=0.754;g=9.81;
L=0.7;mw=13.1;mh=0.86;mk=0.04;m=14.15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chọn n,r, thời gian mô phỏng
tmax=15;
n=9;r=15000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chọn bộ điều khiển + SP 
sp1=1;sp2=0.5;
flag =5;
IS =1;
w_change=0;
l_change=0 ;
e=zeros(1,r);
eh=zeros(1,r);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Chọn thông số cho bộ quan sát trạng thái 
sp3=2;
T_set=2;
b0=1/(13.1);
b0_n=1/(0.86);
T_sam=0.001;
T_sam_n=0.001;
s_cl=-6/T_set;
s_eso=9*s_cl;
z_eso=exp(s_eso*T_sam);
Kd=-2*s_cl;                                
Kp=(s_cl)^2;
% l1=-3*s_eso;
% l2=3*(s_eso)^2;
% l3=-(s_eso)^3;
l1= 1-(z_eso)^3;
l2=3/(2*T_sam)*(1-z_eso)^2*(1+z_eso);
l3=1/T_sam^2*(1-z_eso)^3;
%Bộ quan sát trạng thái
A_LC= [-l1 1 0; -l2 0 1; -l3 0 0];
B=[0 b0 0];
L=[l1;l2;l3];
x1=zeros(1,r);
x2=zeros(1,r);
x3=zeros(1,r);
x1_n=zeros(1,r);
x2_n=zeros(1,r);
x3_n=zeros(1,r);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=tmax/(r-1);
ds=L/(n-1); 
F1=zeros(1,r);
h=zeros(1,r);h(1:3)=ds;
w3=zeros(1,r);
F2=zeros(1,r);
l=2; %% vị trí ban đầu của xe nâng
w=zeros(n,r);

mk=zeros(1,r);mk(:)=0.04;

L=zeros(1,r);L(:)=0.7;

ds=zeros(1,r);ds(:)=0.7/(n-1);
if w_change == 1 
mk(1:round(r/3))=0.04;
mk(round(r/3):round(2*r/3))=2;
mk(round(2*r/3):r)=5;
end 
if l_change == 1 
L= zeros(1,r);
L(1:round(r/3))=0.7;
ds(1:round(r/3))=0.7/n;
L(round(r/3):round(2*r/3))=0.8;
ds(round(r/3):round(2*r/3))=0.8/n;
L(round(2*r/3):r)=1;
ds(round(2*r/3):r)=1/n;
end
%IS
f=6.13415;
 tis1 = 2*3.14/(3*2*3.14*f);
    tis2 = 4*3.14/(3*2*3.14*f);
    tis3 = 2*3.14/(2*3.14*f);
    ris1 = tis1/dt;
    ris2 = tis2/dt;
    ris3 = tis3/dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       MÔ PHỎNG
for j =3:(r-2)
wsss0=(w(3,j)-2*w(2,j)+w(1,j))/(2*ds(j)^3);
w(1,j+1)=(F1(j+1)-EI*wsss0)*(dt^2/mw)+2*w(1,j)-w(1,j-1);%5b
for i=3:(n-2)
S2=(-EI*dt^2)/(ds(j)^4*pA);
S21=-EI*dt^2/ds(j)^3*mh;
wssss=w(i+2,j)-4*w(i+1,j)+6*w(i,j)-4*w(i-1,j)+w(i-2,j);
ddtx3=(w3(j)-2*w3(j-1)+w3(j-2))/(dt^2);
dyx2=(w(l+1,j-1)-w(l-1,j-1))/2*ds(j);
h(j+1)=(F2(j+1)-mh*g-mh*ddtx3*dyx2)*dt^2/mh+2*h(j)-h(j-1);%5c
if h(j+1)<ds(j)
    l=2;
    h(j+1)=ds(j);
end
if h(j+1)>L(j)-ds(j)
    l=n-1;
    h(j+1)=L(j)-ds(j);
end
if h(j+1)> ds(j) && h(j+1)<L(j+1)-ds(j)
    l= ceil(h(j+1)/ds(j));
end
if i == l
w(i,j+1) =S21*wssss+2*w(i,j)-w(i,j-1);%5d
w3(j+1)=w(i,j+1)-w(1,j+1);
end
if i ~= l 
w(i,j+1)=S2*wssss+2*w(i,j)-w(i,j-1);%5a
w3(j+1)=w(i,j+1)-w(1,j+1);
end
end
 S3=(EI*(dt^2))/(mk(j)*2*ds(j)^3);
 w(2,j+1)=w(1,j+1);
 wsssl=(-2*w(n,j)+3*w(n-1,j)-w(n-2,j));
 w(n,j+1)=2*w(n,j)-w(n,j-1)+S3*wsssl;%5e
 w(n-1,j+1)=(w(n,j+1)+w(n-2,j+1))/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flag =0 = check 
% flag = 1 = robust -PD
% flag = 2 = barrier
if flag == 0
F1(1:r/2)=10;
F1(r/2:r)=0;
F2(1:r/2)=10;
check2= zeros(1,r);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                              
if flag == 1 
%ROBUST-PD
%thông số robust - PD
k1=-5;kx=10;
kn=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
dx1=(w(1,j+1)-w(1,j))/dt;
F1(j+2)=k1*(w(1,j+1)-sp1)-kx*dx1;
% if F1(j+2)>50
%     F1(j+2)=50;
% end
% if F1(j+2)<-50
%     F1(j+2)=-50;
% end
S2=(w(l+1,j+1)-2*w(l,j+1)+w(l-1,j+1))/(ds(j)^3);
dx2=(w(l+1,j+1)-w(l-1,j+1))/2*ds(j);
dtx2=(h(j+1)-h(j))/dt;
F2(j+2)=mh*g+EI*S2*dx2-(h(j+1)-sp2)-kn*dtx2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                              
if flag == 2    
% barrier
% thông số 
k1=10;
k2=0;
k3=20;
k0=1;
kc=1;
z=w(l,j+1);
kmax=0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dw1=(w(1,j+1)-w(1,j))/dt;
S1=kc/((kmax^2-z^2)*2.3)+k2;
F1(j+2)=-k1*(w(1,j+1)-sp1)-k3*dw1-S1*kmax*k0*dw1-S1*kmax;
end
if flag == 3
% PD 
% thông số
ki1 =8; kd1 =10 ;
ki2 =1200; kd2 =1200;
%e(j+2)=(sp1-w(1,j+1));
eh(j+1)=(sp2-h(j+1));
dx2=(sp2-h(j+1));
F1(j+2)= ki1*(sp1-w(1,j+1)) + kd1*(w(1,j)-w(1,j+1))/dt;
F2(j+2)= ki2*(sp2-h(j+1)) + kd2*(h(j)-h(j+1))/dt;
end
if flag ==4
% ADRC rời rạc hóa
x1(j+2)=(1-l1)*x1(j+1)+T_sam*(1-l1)*x2(j+1)+T_sam^2/2*(1-l1)*x3(j+1)+(b0*T_sam^2/2-l1*b0*T_sam^2/2)*F1(j+1)+l1*w(1,j+1);
x2(j+2)=-l2*x1(j+1)+(1-l2*T_sam)*x2(j+1)+(T_sam-l2*T_sam^2/2)*x3(j+1)+(b0*T_sam-l2*b0*T_sam^2/2)*F1(j+1)+l2*w(1,j+1);
x3(j+2)=-l3*x1(j+1)-l3*T_sam*x2(j+1)+(1-l3*T_sam^2/2)*x3(j+1)-l3*b0*T_sam^2/2*F1(j+1)+l3*w(1,j+1);
% Lực điều khiển F1 với giới hạn
F1(j+2) = (1/b0)*((Kp*(sp1 - x1(j+2)) - Kd*x2(j+2)) - x3(j+2));
F1(j+2) = min(max(F1(j+2), -50), 50);  % Giới hạn lực F1
%
x1_n(j+2)=(1-l1)*x1_n(j+1)+T_sam_n*(1-l1)*x2_n(j+1)+T_sam_n^2/2*(1-l1)*x3_n(j+1)+(b0_n*T_sam_n^2/2-l1*b0_n*T_sam_n^2/2)*F2(j+1)+l1*h(j+1);
x2_n(j+2)=-l2*x1_n(j+1)+(1-l2*T_sam_n)*x2_n(j+1)+(T_sam_n-l2*T_sam_n^2/2)*x3_n(j+1)+(b0_n*T_sam_n-l2*b0_n*T_sam_n^2/2)*F2(j+1)+l2*h(j+1);
x3_n(j+2)=-l3*x1_n(j+1)-l3*T_sam_n*x2_n(j+1)+(1-l3*T_sam_n^2/2)*x3_n(j+1)-l3*b0_n*T_sam_n^2/2*F2(j+1)+l3*h(j+1);
% Lực điều khiển F2 với giới hạn
F2(j+2) = (1/b0_n)*((Kp*(sp2 - x1_n(j+2)) - Kd*x2_n(j+2)) - x3_n(j+2));
 F2(j+2) = min(max(F2(j+2), -50), 50);  % Giới hạn lực F2
end
if flag ==5
    % ADRC rời rạc hóa + IS
    if IS==1
     spIS1= sp1*(1/2*us(j*dt)+1/2*us(j*dt-0.081511));
    end
     if IS==2
     spIS1= sp1*(1/4*us(j*dt)+1/2*us(j*dt-0.081511)+1/4*us(j*dt-0.163022));
     end
     if IS ==3
     spIS1= sp1*(1/6*us(j*dt)+1/3*us(j*dt-0.054341)+1/3*us(j*dt-0.108681)+1/6*us(j*dt-0.163022));
     end
     %
x1(j+2)=(1-l1)*x1(j+1)+T_sam*(1-l1)*x2(j+1)+T_sam^2/2*(1-l1)*x3(j+1)+(b0*T_sam^2/2-l1*b0*T_sam^2/2)*F1(j+1)+l1*w(1,j+1);
x2(j+2)=-l2*x1(j+1)+(1-l2*T_sam)*x2(j+1)+(T_sam-l2*T_sam^2/2)*x3(j+1)+(b0*T_sam-l2*b0*T_sam^2/2)*F1(j+1)+l2*w(1,j+1);
x3(j+2)=-l3*x1(j+1)-l3*T_sam*x2(j+1)+(1-l3*T_sam^2/2)*x3(j+1)-l3*b0*T_sam^2/2*F1(j+1)+l3*w(1,j+1);
% Lực điều khiển F1 với giới hạn
F1(j+2) = (1/b0)*((Kp*(spIS1 - x1(j+2)) - Kd*x2(j+2)) - x3(j+2));
% F1(j+2) = min(max(F1(j+2), -50), 50);  % Giới hạn lực F1
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
subplot(2,2,4)
grid on
plot(x,w3,LineWidth=1)
title('Độ lắc xe nâng');
xlabel('Thời gian(s)');
ylabel('m');

figure(2)
subplot(1,2,1)
plot(x,F1(1:r),LineWidth=1);
title('F1')
grid on
subplot(1,2,2)
plot(x,F2(1:r),LineWidth=1);
title('F2');
grid on
% figure(3)
% hold on
% grid on
% plot(x,mk,Linewidth=1);
% title("Khối lượng xe")
% figure(4)
% plot(x,L,Linewidth=1);
% title("Chiều dài của thang")
%%Phân tích phổ
figure(4)
Fs=1/dt;
time=0:dt:(r-2)*dt-dt;
l1=r;
fft_w=fft(w(n,:)-w(1,:),l1)*(2/l1);
abs_w=abs(fft_w);
freq=0:(1/time(end)):Fs/2-(1/time(end));
plot(freq,abs_w(1:length(freq)),'LineWidth',0.8);
hold on
grid on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 9);
xlabel('Frequency(Hz)');
ylabel('Amplitude(m)');
title('Pho dao dong');
