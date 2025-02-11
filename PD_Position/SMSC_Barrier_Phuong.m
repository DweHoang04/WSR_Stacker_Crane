clc; clear;

pA=0.297; EI=0.754; g=9.81;
L=0.7; mw=13.1; mh=0.86; mk=0.04; m=14.15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chọn n,r, thời gian mô phỏng
tmax=15;
n=9;r=10000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chọn bộ điều khiển + SP 
sp1=1; sp2=0.5;
flag = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=tmax/(r-1);
ds=L/(n-1); 
F1=zeros(1,r);
h=zeros(1,r); h(1:3)=ds;
w3=zeros(1,r);
F2=zeros(1,r);
l=2; 
w=zeros(n,r);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       MÔ PHỎNG
for j =3:r-1
    wsss0=(w(3,j)-2*w(2,j)+w(1,j))/(2*ds^3);
    w(1,j+1)=(F1(j+1)-EI*wsss0)*(dt^2/mw)+2*w(1,j)-w(1,j-1);
    for i=3:n-2
        S2=(-EI*dt^2)/(ds^4*pA);
        S21=-EI*dt^2/ds^3*mh;
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
        if i == l
        w(i,j+1) =S21*wssss+2*w(i,j)-w(i,j-1);
        w3(j+1)=w(i,j+1)-w(1,j+1);
        end
        if i ~= l 
        w(i,j+1)=S2*wssss+2*w(i,j)-w(i,j-1);
        w3(j+1)=w(i,j+1)-w(1,j+1);
        end
     end
     S3=(EI*(dt^2))/(mk*2*ds^3);
     w(2,j+1)=w(1,j+1);
     wsssl=(-2*w(n,j)+3*w(n-1,j)-w(n-2,j));
     w(n,j+1)=2*w(n,j)-w(n,j-1)+S3*wsssl;
     w(n-1,j+1)=(w(n,j+1)+w(n-2,j+1))/2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                              
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % flag = 0 = check 
    % flag = 1 = robust - PD
    % flag = 2 = barrier
    if flag == 0
        F1(1:r/2)=10;
        F1(r/2:r)=0;
        F2(1:r/2)=0;
        check2= zeros(1,r);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                              
    if flag == 1 
        %ROBUST-PD
        %thông số robust - PD
        k1 = -2; kx = 10;
        kn = 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        dx1=(w(1,j+1)-w(1,j))/dt;
        F1(j+2)=k1*(w(1,j+1)-sp1)-kx*dx1;
        if F1(j+2)>50
            F1(j+2)=50;
        end
        if F1(j+2)<-50
            F1(j+2)=-50;
        end
        S2=(w(l+1,j+1)-2*w(l,j+1)+w(l-1,j+1))/(ds^3);
        dx2=(w(l+1,j+1)-w(l-1,j+1))/(2*ds);
        dtx2=(h(j+1)-h(j))/dt;
        F2(j+2)=mh*g+EI*S2*dx2-(h(j+1)-sp2)-kn*dtx2;
        if F2(j+2)>50
            F2(j+2)=50;
        end
        if F2(j+2)<-50
            F2(j+2)=-50;
        end
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
plot(x,w(n,:)-w(1,:),'b',LineWidth=1);
title('Độ lắc điểm cuối');
xlabel('Thời gian(s)');
ylabel('m');
check = w3-w(1,:);

subplot(2,2,3)
hold on
grid on
plot(x,h,'r');
title('Vị trí xe nâng',LineWidth=1);
xlabel('Thời gian(s)');
ylabel('m');
title('Vị trí xe nâng');

subplot(2,2,4)
grid on
hold on
plot(x,w3,LineWidth=1)
title('Độ lắc xe nâng');
xlabel('Thời gian(s)');
ylabel('m');

figure(2)
subplot(2,1,1);
grid on;
hold on;
plot(x,F1(1:r));
title({'Lực F1'});
ylabel('Lực F1 (N)','FontSize',12);
xlabel('Thời gian (s)','FontSize',12);