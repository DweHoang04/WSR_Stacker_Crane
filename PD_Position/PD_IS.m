clc; close all
% 3 thong so can chinh: k1, k2, k3, (k0 có thể chỉnh nếu muốn)
% Các tham số tăng:
% k1: giảm thời gian xác lập, tăng quá độ điều chỉnh
% k2: giảm rung nhanh hơn, loe nhanh hơn
% k3: loe nhanh hơn, vị trí xe đáp ứng chậm hơn, không quá độ
% k0: có thể chỉnh nhưng chưa rõ hướng thay đổi
% c: để giữ tại 0.1 để khớp thực nghiệm

set(0, 'DefaultAxesFontName', 'Times New Roman'); % Phông chữ cho trục
set(0, 'DefaultTextFontName', 'Times New Roman'); % Phông chữ cho văn bản
set(0, 'DefaultAxesFontSize', 10); % Kích thước chữ cho trục
set(0, 'DefaultTextFontSize', 10); % Kích thước chữ cho văn bản
k1=250;k3=130;xd=0.2; SP=xd;

pA=2.4882;EI=1.05;
l=0.6;
mw=13.1;
mk=0.2;
tmax=15;
r=140000;
n = 40;
dt=tmax/(r-1);
heso = 1;
c=0.11;

IS = 0; %% (0: khong dung IS; 1: ZV; 2: ZVD; 3: ETM4)
%%Tinh toan cho IS
f = 1.73336;

if IS == 1 %ZV
    tis1 = 3.14/(2*3.14*f);
    ris1 = tis1/dt;
elseif IS == 2 %ZVD
    tis1 = 3.14/(2*3.14*f);
    tis2 = 2*3.14/(2*3.14*f);
    ris1 = tis1/dt;
    ris2 = tis2/dt;
elseif IS == 3 %ETM4
    tis1 = 2*3.14/(3*2*3.14*f);
    tis2 = 4*3.14/(3*2*3.14*f);
    tis3 = 2*3.14/(2*3.14*f);
    ris1 = tis1/dt;
    ris2 = tis2/dt;
    ris3 = tis3/dt;
end
ds=l/(n-1);
S1=(EI*(dt^2))/(mw*(ds^3));
S2=(-EI*(dt^2))/((ds^4)*pA);
S3=(EI*(dt^2))/(mk*ds^3);
%nhap input
%input gom 2 phan = F1(t)+F2(t)
F1=zeros(1,r);
% F1(1:50000)=1;
% F1(50001:100000)=-1;
% F1(100001:500000)=0;
w=zeros(n,r);
A = [-0.5 30 -0.1e-5 2e-2];
wt=zeros(n,r);
wsss0=zeros(1,r);
wssst0=zeros(1,r);
wt0=zeros(1, r);
wl=zeros(1,r);
wtl=zeros(1,r);
SP212=zeros(1,r);
SP212(r-1) = xd;
SP212(r) = xd;
%for j =1:100000
for j =1:r-2
    %     t = j*dt;
    %     xd = QD212(1,0.5,1,SP,t);
    if IS == 1
        if j<ris1
            xd=1/2*SP;
        else 
            xd=1*SP;
        end
    elseif IS == 2
        if j<ris1
            xd=1/4*SP;
        elseif j<ris2
            xd=(1/4+1/2)*SP;
        else
            xd=1*SP;
        end
    elseif IS == 3
        if j<ris1
            xd=1/6*SP;
        elseif j<ris2
            xd=(1/6+1/3)*SP;
        elseif j<ris3
            xd=(1/6+1/3+1/3)*SP;
        else
            xd=1*SP;
        end
    end

    SP212(j)= xd;
    wsss0(j+1)=(w(3,j+1)-2*w(2,j+1)+w(1,j+1))/(ds^3);
    wssst0(j+1)=(wsss0(j+1)-wsss0(j))/dt;
    wt0(j+1)=(w(1,j+1)-w(1,j))/dt;
    wtl(j+1)=(wl(j+1)-wl(j))/dt;
    F1(j+2)=-(k1)*(w(1,j+1)-xd)-(k3)*(w(1,j+1)-w(1,j))/dt;
    w(1,j+2)=(F1(j+2)-EI*wsss0(j+2))*(dt^2/mw)+2*w(1,j+1)-w(1,j);
    % w(2,j+2)=(w(1,j+2)+w(3,j+2))/2;
    for i=3:n-2
        S2=(-EI*dt^2)/(ds^4*pA);

        wt(n,j+1)=(w(n,j+1)-w(n,j))/dt;
        wssss=w(i+2,j+1)-4*w(i+1,j+1)+6*w(i,j+1)-4*w(i-1,j+1)+w(i-2,j+1);
        w(i,j+2)=S2*wssss+2*w(i,j+1)-w(i,j);
        w(i,j+2)=(S2*wssss+(2+c*dt/pA)*w(i,j+1)-w(i,j))/(1+c*dt/pA);

    end
    w(2,j+2)=(w(1,j+2)+w(3,j+2))/2;
    wsssl=(-2*w(n,j+1)+3*w(n-1,j+1)-w(n-2,j+1));
    w(n,j+2)=2*w(n,j+1)-w(n,j)+S3*wsssl;
    w(n-1,j+2)=(w(n,j+2)+w(n-2,j+2))/2;
    wl(j+2)=w(n,j+2)-w(1,j+2);
end
% w=w*1000;
x=0:dt:tmax;
figure('Position', [100, 100, 600, 300])
figure(1)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 9);
hold on
grid on
plot(x, F1,'color', [0.07,0.62,1.00],'linewidth',2);
if IS==0
    F10 = F1;
elseif IS ==1
    F11 = F1;
elseif IS ==2
    F12 = F1;
else
    F13 = F1;
end
xlabel('Thời gian (s)', 'FontSize',11);
ylabel('Lực F_1 (N)', 'FontSize',11);
axis([0 10 -10 55]);

%title('F1')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n=5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             %FIGURE %
x=0:dt:tmax;
figure('Position', [100, 100, 600 300])

figure(2)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 9);
grid on
hold on
plot(x,SP212,'--r','linewidth',2);
plot(x,w(1,:),'color', [0.07,0.62,1.00],'linewidth',2);
if IS==0
    x10 = w(1,:);
elseif IS ==1
    x11 = w(1,:);
elseif IS ==2
    x12 = w(1,:);
else
    x13 = w(1,:);
end


%plot(x,w1(1,:), '--', LineWidth=0.7);
%legend('Anti-vibration position control','Position control');
%title('Vị trí xe con')
legend('SP','PV');
xlabel('Thời gian (s)', 'FontSize',11);
ylabel('Vị trí xe con (m)', 'FontSize',11);
axis([0 10 0 0.25]);

% xticks([0 2 4 6 8 10]);
% yticks([0 0.05 0.1 0.15 0.2 0.25]);
%
%
% figure(2)%diem giua
% hold on
% grid on
% plot(x,w(ceil(n/2),:)-w(1,:),'b');
% plot(x,w1(ceil(n1/2),:)-w1(1,:),'r');
% plot(x,w3(ceil(n3/2),:)-w3(1,:),'y');
% plot(x,w4(ceil(n4/2),:)-w4(1,:),'g');
% legend(num2str(n),num2str(n1), num2str(n3), num2str(n4));
% title('Diem giua');

figure('Position', [100, 100, 600 300])
figure(3)%diem cuoi
set(gca, 'FontName', 'Times New Roman', 'FontSize', 9);
hold on
grid on
plot(x,(w(n,:)-w(1,:))*1000, 'color', [0.07,0.62,1.00],'linewidth',2);
% plot(x,(w(n,:)-w(1,:))*1000, 'color', [0.47,0.62,1.00],'linewidth',2, 'r');

% plot(x,(w(n,:)-w(1,:))*1000, 'r','linewidth',2);

if IS==0
    d10 = w(n,:)-w(1,:);
elseif IS ==1
    d11 = w(n,:)-w(1,:);
elseif IS ==2
    d12 = w(n,:)-w(1,:);
else
    d13 = w(n,:)-w(1,:);
end
xlabel('Thời gian (s)', 'FontSize',11);
ylabel('Dao động (mm)', 'FontSize',11);
axis([0 10 -40 40]);


% xticks([0 2 4 6 8 10]);
% yticks([-0.08 -0.04 0 0.04 0.08]);

% %Phan tich pho%%
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



n1 = 5; r1 = 3500;
wg = zeros(n1,r1);
rn1 = n/n1; rr1 = r/r1;
for i=1:n1
    for j=1:r1
        wg(i,j) = w(rn1*i-(rn1-1), rr1*j-(rr1-1));
    end
end
xc = (wg-wg(1,:))*1000;
xc1 = xc';

figure('Position', [700 100 600 300]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 9);
[X, Y] = meshgrid(0:(tmax/(r1-1)):tmax,0:(l/(n1-1)):l);
meshc(Y,X,xc);
% ylabel('t(s)');
% xlabel('Y(m)');
% zlabel('ω(Y,t)(mm)');

ylabel('$t$(s)', 'Interpreter', 'latex', 'FontSize',9);
xlabel('$Y$(m)', 'Interpreter', 'latex', 'FontSize',9);
zlabel('$\omega(Y,t)$(mm)', 'Interpreter', 'latex', 'FontSize',9);
view(65,10);
axis([0 0.6 0 10 -40 40])


figure(Position=[200 200 600 300])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 9);
subplot(2,2,1)
% plot(x,w(1,:),'color', [0.07,0.62,1.00],'linewidth',2);
plot(x,w(1,:),'color', [0, 0, 0],'linewidth',1);

grid on
hold on
xlabel('$t$(s)', 'Interpreter', 'latex', 'FontSize',9);
ylabel('$x_1$(m)','Interpreter', 'latex', 'FontSize',9);
axis([0 10 0 0.25]);
subplot(2,2,2)
% plot(x,(w(n,:)-w(1,:))*1000, 'color', [0.07,0.62,1.00],'linewidth',2);
plot(x,(w(n,:)-w(1,:))*1000, 'color', [0 0 0],'linewidth',1);
xlabel('$t$(s)', 'Interpreter', 'latex', 'FontSize',9);
ylabel('$\omega(L,t)$(mm)', 'Interpreter', 'latex', 'FontSize',9);

axis([0 10 -40 40]);
hold on; grid on;
subplot(2,2,3)
% plot(x, F1,'color', [0.07,0.62,1.00],'linewidth',2);
plot(x, F1,'color', [0 0 0],'linewidth',1);
xlabel('$t$(s)', 'Interpreter', 'latex');
ylabel('$F_1$(N)', 'Interpreter', 'latex');
axis([0 10 -10 55]);
hold on
grid on
subplot(2,2,4)
plot(freq,1000*abs_w(1:length(freq)),'k','linewidth',1);
xlabel('$f$(Hz)', 'Interpreter', 'latex');
ylabel('$A$(mm)', 'Interpreter', 'latex');
% ylabel('Am(mm)', 'FontSize',9);
axis([0 8 0 40]);
grid on; hold on;

figure()
subplot(2,2,1)
% plot(x,w(1,:),'color', [0.07,0.62,1.00],'linewidth',2);
plot(x,w(1,:),'color', [0, 0, 0],'linewidth',1, 'LineStyle','--');
legend('IS+PD', 'PD');

subplot(2,2,2)
% plot(x,(w(n,:)-w(1,:))*1000, 'color', [0.07,0.62,1.00],'linewidth',2);
plot(x,(w(n,:)-w(1,:))*1000, 'color', [0 0 0],'linewidth',1,  'LineStyle','--');
legend('IS+PD', 'PD');

axis([0 10 -40 40]);
hold on; grid on;
subplot(2,2,3)
% plot(x, F1,'color', [0.07,0.62,1.00],'linewidth',2);
plot(x, F1,'color', [0 0 0],'linewidth',1, 'LineStyle','--');
legend('IS+PD', 'PD');
axis([0 10 -10 55]);
hold on
grid on

subplot(2,2,4)
plot(x,(w(19,:)-w(1,:))*1000, 'color', [0 0 0],'linewidth',1, 'LineStyle','--');
legend('IS+PD', 'PD');
% ylabel('Am(mm)', 'FontSize',9);
axis([0 8 -40 40]);
grid on; hold on;

figure(6)
subplot(2,2,1)
% plot(x,w(1,:),'color', [0.07,0.62,1.00],'linewidth',2);
plot(x,w(1,:),'color', [0, 0, 0],'linewidth',1, 'LineStyle','--');
legend('PDv1', 'PD');
grid on; hold on;

subplot(2,2,2)
% plot(x,(w(n,:)-w(1,:))*1000, 'color', [0.07,0.62,1.00],'linewidth',2);
plot(x,(w(n,:)-w(1,:))*1000, 'color', [0 0 0],'linewidth',1,  'LineStyle','--');
legend('PDv1', 'PD');

axis([0 10 -40 40]);
hold on; grid on;
subplot(2,2,3)
% plot(x, F1,'color', [0.07,0.62,1.00],'linewidth',2);
plot(x, F1,'color', [0 0 0],'linewidth',1, 'LineStyle','--');
legend('PDv1', 'PD');
axis([0 10 -10 55]);
hold on
grid on

subplot(2,2,4)
plot(x,(w(19,:)-w(1,:))*1000, 'color', [0 0 0],'linewidth',1, 'LineStyle','--');
legend('PDv1', 'PD');
% ylabel('Am(mm)', 'FontSize',9);
axis([0 8 -40 40]);
grid on; hold on;
