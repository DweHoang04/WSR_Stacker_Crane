clc;clear;close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pA=0.297;EI=0.754;g=9.81;
L=0.7;mw=13.1;mh=0.86;mk=0.04;m=14.15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chọn n,r, thời gian mô phỏng
tmax=15;
n=9;r=14000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chọn bộ điều khiển + SP 
% sp1=1;sp2=0.2;
sp1=1;sp2=0.3;
flag = 2;
w_change=0;
l_change=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BIẾN CHO ĐIỀU KHIỂN BOUNDARY
d11=zeros(1,r+1);d21=zeros(1,r+1);d31=zeros(1,r+1);
dF1=zeros(1,r+1);a4=zeros(1,r+1);a1=zeros(1,r+1);
u1=zeros(1,r+1);a2=zeros(1,r+1);
v1=zeros(1,r+1);
U1=zeros(1,r+1);
% BIẾN CHO ĐIỀU KHIỂN INPUT CONS
p1=zeros(1,r);
Fc1=zeros(1,r);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=tmax/(r-1);
ds=L/(n-1); 
e=zeros(1,r);
F1=zeros(1,r);
h=zeros(1,r);h(1:3)=ds;
w3=zeros(1,r);
F2=zeros(1,r);
l=2; %% vị trí ban đầu của xe nâng
w=zeros(n,r);

mk=zeros(1,r);mk(:)=0.04;

L=zeros(1,r);L(:)=0.7;

z1=zeros(1,r);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       MÔ PHỎNG
for j =3:r-1
wsss0=(w(3,j)-2*w(2,j)+w(1,j))/(2*ds(j)^3);
 w(1,j+1)=(F1(j+1)-EI*wsss0)*(dt^2/mw)+2*w(1,j)-w(1,j-1);
 % w(1,j+1)=(F1(j+1)-EI*wsss0+2*sin(j*pi/4))*(dt^2/mw)+2*w(1,j)-w(1,j-1);

for i=3:n-2
S2=(-EI*dt^2)/(ds(j)^4*pA);
S21=(-EI*dt^2)/ds(j)^3*mh;
c=0.11;
C=-c*dt^2*(w(i,j)-w(i,j-1))/dt;
wssss=w(i+2,j)-4*w(i+1,j)+6*w(i,j)-4*w(i-1,j)+w(i-2,j);
ddtx3=(w3(j)-2*w3(j-1)+w3(j-2))/(dt^2);
dyx2=(w(l+1,j-1)-w(l-1,j-1))/2*ds(j);
h(j+1)=(F2(j+1)-mh*g-mh*ddtx3*dyx2)*dt^2/mh+2*h(j)-h(j-1);
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
 S3=(EI*(dt^2))/(mk(j)*2*ds(j)^3);
 w(2,j+1)=w(1,j+1);
 wsssl=(-2*w(n,j)+3*w(n-1,j)-w(n-2,j));
 w(n,j+1)=2*w(n,j)-w(n,j-1)+S3*wsssl;
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
k1=10;kx=20;kn=1.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% e(j+1)=(w(1,j+1)-w(1,j));
% de=(e(j+1)-e(j))/dt;
dx1=(w(1,j+1)-w(1,j))/dt;
F1(j+2)=-k1*(w(1,j+1)-sp1)-kx*dx1;
% if F1(j+2)>50
%     F1(j+2)=50;
% end
% if F1(j+2)<-5
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
k1=20;kmax=0.004;k2=0.01;k3=30;k0=10;kc=0.001;kn=1.5;
z=w(l,j+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dwl=(w(1,j+1)-w(1,j))/dt;
S1=kc/((kmax^2-z^2)*2.3)+k2;
F1(j+2)=-k1*(w(1,j+1)-sp1)-k3*dwl-S1*kmax*k0*dwl-S1*kmax;

S2=(w(l+1,j+1)-2*w(l,j+1)+w(l-1,j+1))/(ds(j)^3);
dx2=(w(l+1,j+1)-w(l-1,j+1))/2*ds(j);
dtx2=(h(j+1)-h(j))/dt;
F2(j+2)=mh*g+EI*S2*dx2-(h(j+1)-sp2)-kn*dtx2;

end
if flag == 3
% PD 
% thông số
ki1 = 1; kd1 = 10 ;
ki2 = 0; kd2 = 0;
e(j+1)=(sp1-w(1,j+1));
de=(e(j+1)-e(j))/dt;
F1(j+2)= ki1*(sp1-w(1,j+1)) + kd1*de;


end
if flag ==4
 % thong so 
 a1=1;k1=25;k2=15;kn=1.5;

 % disturbance observer
dwl=(w(1,j+1)-w(1,j))/dt;
dz1=-a1*(F1(j+1)-EI*wsss0)-a1*(z1(j+1)+a1*mw*dwl);
z1(j+2)=dz1*dt+z1(j+1);
d1e = z1(j+2)+a1*mw*dwl;
F1(j+2)=-k1*dwl-d1e-k2*(w(1,j+1)-sp1);

S2=(w(l+1,j+1)-2*w(l,j+1)+w(l-1,j+1))/(ds(j)^3);
dx2=(w(l+1,j+1)-w(l-1,j+1))/2*ds(j);
dtx2=(h(j+1)-h(j))/dt;
F2(j+2)=mh*g+EI*S2*dx2-(h(j+1)-sp2)-kn*dtx2;

end
% if flag == 5
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% um=1;vm=25;beta=1;
% c1=1;c7=4;c4=1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dwl=(w(n,j+1)-w(n,j))/dt;
% d11(j+1)=dwl+beta*L(j+1)*(w(n,j+1)-w(n-1,j+1))/ds(j+1);
% dsdt=(w(n,j+1)-w(n-1,j+1)-w(n,j)+w(n-1,j))/(ds(j+1)*dt);
% a1(j+1)=-mw*beta*L(j+1)*dsdt-2*c1*dwl;
% da1=(a1(j+1)-a1(j))/dt;
% d21(j+1)=F1(j+1)-a1(j+1);
% a4(j+1)=-1/mw*d11(j+1)-1/mw*c4*d21(j+1)*da1;
% d31(j+1)=dF1(j+1)-a4(j+1);
% da4=(a4(j+1)-a4(j))/dt;
% U1(j+2)=-1/mw*c7*d31(j+1)-d21(j+1)+da4;
% dv1= (1/(cosh(v1(j+1)/vm))^2)*U1(j+2);
% v1(j+2)=dt*dv1+v1(j+1);
% dF1(j+2)= vm*tanh(v1(j+2)/vm);
% % F1(j+2)=dF1(j+2)*dt+F1(j+1);
% du1=(1/(cosh(u1(j+1)/vm))^2)*dF1(j+2);
% u1(j+2)=du1*dt+u1(j+1);
% F1(j+2)=um*tanh(u1(j+2)/um);
% end
if flag ==5 
um=15;vm=60;
c1=10;c3=13;c2=40;k1=0.12;kn=1.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dw1=(w(1,j+1)-w(1,j))/dt;
wsss=(w(3,j+1)-3*w(2,j+1)+2*w(1,j+1))/ds(j+1)^3;

e(j+1)=(sp1-w(1,j+1));
de=(e(j+1)-e(j))/dt;

a1(j+1)=-c1*dw1+1/2*EI*wsss+1/2*k1*e(j+1);
d11=dw1;
d21=F1(j+1)-a1(j+1);
a2(j+1)=-2*(1/mw)*d11-(1/mw)*c2*d21+(a1(j+1)-a1(j))/dt;
d31(j+1)=dF1(j+1)-a2(j+1);
da2=(a2(j+1)-a2(j))/dt;
U1(j+2)=(-1/mw)*c3*d31(j+1)-d21+da2;
dv1= (1/(cosh(v1(j+1)/vm))^2)*U1(j+2);
v1(j+2)=dt*dv1+v1(j+1);
dF1(j+2)= vm*tanh(v1(j+2)/vm);
du1=(1/(cosh(u1(j+1)/vm))^2)*dF1(j+2);
u1(j+2)=du1*dt+u1(j+1);
F1(j+2)=um*tanh(u1(j+2)/um);

S2=(w(l+1,j+1)-2*w(l,j+1)+w(l-1,j+1))/(ds(j)^3);
dx2=(w(l+1,j+1)-w(l-1,j+1))/2*ds(j);
dtx2=(h(j+1)-h(j))/dt;
F2(j+2)=mh*g+EI*S2*dx2-(h(j+1)-sp2)-kn*dtx2;
end

if flag == 6 
    
    %% thong so 
    kmax=0.01;
    g1=0.1;
    k1=1;
    c1=1;
    k11=1;
    k2=1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        wsss0=(w(3,j)-3*w(2,j)+2*w(1,j));
        wsssl=(-2*w(n,j)+3*w(n-1,j)-w(n-2,j));
        dw1=(w(1,j+1)-w(1,j))/dt;
        z(j+1)=w(n,j+1)-w(1,j+1);
        a = log((2 * kmax^2) / (kmax^2 - z(j+1)^2))/log(2.71);
        dz=(z(j+1)-z(j))/dt;
        s1=EI*wsss0/ds(j)^3;
      v1(j+2)=-k1*dz-s1-k11*dz/a-mw*dz/a*z(j+1)*dz/(kmax^2-z(j+1)^2)-(k2*(w(1,j)-sp1))/(a)-mw*EI*wsssl/(ds(j)^3*mk(j)*a)  ;
       %v1(j+2)=-k1*dz-s1-k11*dz/a-mw*dz/a*z(j+1)*dz/(kmax^2-z(j+1)^2);
        %%%% uoc luong sai so
        p1(j+2)=(g1*dz*v1(j+2)*a-p1(j+1))*dt+p1(j+1);
        Fc1(j+2)=p1(j+2)*v1(j+2);
        F1(j+2)=c1*Fc1(j+2);
end
end
end
x=0:dt:tmax;
y=0:ds:L;
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
for j=1:r
dolac(j)=w(n,j)-w(1,j);
end
plot(x,dolac,'b',LineWidth=1);
title('Độ lắc điểm cuối');
xlabel('Thời gian(s)');
ylabel('m');
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
check=dolac-w3;
 subplot(1,2,2)
plot(x,F2(1:r),LineWidth=1);
title('F2');
grid on

n1 = 3; r1 = 1000;
w = zeros(n,r);
rn1 = n/n1; rr1 = r/r1;
% for i=1:n1
%     for j=1:r1
%         wg(i,j) = w(rn1*i-(rn1-1), rr1*j-(rr1-1));
%     end
% end
xc = (w-w(1,:))*1000;
xc1 = xc';

figure('Position', [700 100 600 300]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 9);
[X, Y] = meshgrid(0:(tmax/(r-1)):tmax,0:(l/(n-1)):l);
meshc(Y,X,xc);
% ylabel('t(s)');
% xlabel('Y(m)');
% zlabel('ω(Y,t)(mm)');

ylabel('$t$(s)', 'Interpreter', 'latex', 'FontSize',9);
xlabel('$Y$(m)', 'Interpreter', 'latex', 'FontSize',9);
zlabel('$\omega(Y,t)$(mm)', 'Interpreter', 'latex', 'FontSize',9);
view(65,10);
axis([0 0.6 0 10 -40 40])

% figure(3)
 
% hold on
% grid on
% plot(x,mk,Linewidth=1);
% title("Khối lượng xe")
% figure(4)
% plot(x,L,Linewidth=1);
% title("Chiều dài của thang")

