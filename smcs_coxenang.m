clc; clear;
pA=0.297;EI=0.754;g=9.81;
L=0.63;
mw=13.1;
mh=0.86;
mk=0.04;
tmax=15;
n=5;
r=5000;
dt=tmax/(r-1);
ds=L/(n-1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F1=zeros(1,r);
h=zeros(1,r);
h(1:2)=ds;
w3=zeros(1,r);
F2=zeros(1,r);
l=2;
w=zeros(n,r);
F1(1:r/2)=10;
F2(1:r/2)=mh*g-0.2;
check2= zeros(1,r);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j =3:r-1
wsss0=(w(3,j)-2*w(2,j)+w(1,j))/(2*ds^3);
w(1,j+1)=(F1(j+1)-EI*wsss0)*(dt^2/mw)+2*w(1,j)-w(1,j-1); % 5b
for i=3:n-2
S2=(-EI*dt^2)/(ds^4*pA);
S21=-EI*dt^2/ds^3*mh;
wssss=w(i+2,j)-4*w(i+1,j)+6*w(i,j)-4*w(i-1,j)+w(i-2,j);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ddtx3=(w3(j)-2*w3(j-1)+w3(j-2))/(dt^2);
dyx2=(w(l+1,j-1)-w(l-1,j-1))/2*ds;
check2(j)=ddtx3*dyx2;
h(j+1)=(F2(j+1)-mh*g-mh*ddtx3*dyx2)*dt^2/mh+2*h(j)-h(j-1);
if h(j+1)<ds
    l=2;
    h(j+1)=ds;
elseif h(j+1)>L-ds
    l=n-1;
    h(j+1)=L-ds;
else 
    l= ceil(h(j+1)/ds);
end
if i == l
w(i,j+1) = S21*wssss+2*w(i,j)-w(i,j-1);
w3(j+1)=w(i,j+1)-w(1,j+1);
else
w(i,j+1)=S2*wssss+2*w(i,j)-w(i,j-1);
w3(j+1)=w(i,j+1)-w(1,j+1);
end
end
 S3=(EI*(dt^2))/(mk*2*ds^3);
 w(2,j+1)=w(1,j+1);
 wsssl=(-2*w(n,j)+3*w(n-1,j)-w(n-2,j));
 w(n,j+1)=2*w(n,j)-w(n,j-1)+S3*wsssl; % 5e
 w(n-1,j+1)=(w(n,j+1)+w(n-2,j+1))/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=0:dt:tmax;
figure(1)
grid on 
hold on
plot(x,w(1,:),'b',LineWidth=2.5);
xlabel('Thời gian(s)')
ylim([0 2])
ylabel('Vị trí xe con(m)')
title('Vị trí xe con')

figure(2)
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
check = w3-w(1,:);

figure(3)
plot(x,h,'r');
title('Vị trí xe nâng',LineWidth=1);
xlabel('Thời gian(s)');
ylabel('m');
title('Vị trí xe nâng');
figure(4)
plot(x,w3,LineWidth=1)
title('Độ lắc xe nâng');
xlabel('Thời gian(s)');
ylabel('m');



