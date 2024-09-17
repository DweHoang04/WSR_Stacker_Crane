clc;clear;
k=1;
pA=2.4;EI=0.4;g=9.81;
L=0.7;
mw=13.1;
mh=0.86;
mk=0.2;
m=14.15;
tmax=15;
n=9;
r=5000;
dt=tmax/(r-1);
g=9.8;
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
F1(r/2:r)=0;
F2(1:1000)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j =2:r-1
wsss0=(w(3,j)-2*w(2,j)+w(1,j))/(2*ds^3);
w(1,j+1)=(F1(j+1)-EI*wsss0)*(dt^2/mw)+2*w(1,j)-w(1,j-1);
 for i=3:n-2
S2=(-EI*dt^2)/(ds^4*pA);
wssss=w(i+2,j)-4*w(i+1,j)+6*w(i,j)-4*w(i-1,j)+w(i-2,j);
w(i,j+1)=S2*wssss+2*w(i,j)-w(i,j-1);
 end
 S3=(EI*(dt^2))/(mk*2*ds^3);
 w(2,j+1)=w(1,j+1);
 wsssl=(-2*w(n,j)+3*w(n-1,j)-w(n-2,j));
 w(n,j+1)=2*w(n,j)-w(n,j-1)+S3*wsssl;
 w(n-1,j+1)=(w(n,j+1)+w(n-2,j+1))/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ddtx3=(w3(j+1)-2*w3(j)+w3(j-1))/(dt^2);
dyx2=(w(l,j)-w(l-1,j))/ds;
h(j+1)=(F2(j+1)-mh*g-mh*ddtx3*dyx2)*dt^2/mh+2*h(j)-h(j-1);
l=ceil(h(j+1)/ds);
if h(j+1)<ds
    l=2;
    h(j+1)=ds;
end
if h(j+1)>L-ds
    l=n-1;
    h(j+1)=L-ds;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if l == 2 
    w3(j+1)=w(1,j+1);
end
if l > 2 && l < n -1
 w3(j+1)=-(2*EI*dt^2)*(w(l+1,j+1)-2*w(l,j+1)+w(l-1,j+1))/(ds^3*mh)+2*w3(j)-w3(j-1);
end
if l==n
 w3(j+1)=w(n-1,j+1);
end
check = w3-w(1,:);
end
x=0:dt:tmax;
figure(1)
grid on 
hold on
plot(x,w(1,:),'b',LineWidth=2.5);
xlabel('Thời gian(s)')
ylim([0 2])
ylabel('Vị trí xe con(m)')
title('Vị trí xe con')
figure(3)
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








