clc; clear; close all;
% Thông số
pA = 0.297; EI = 0.754; g = 9.81;
L = 0.7; mw = 13.1; mh = 0.86; mk = 0.04; m = 14.15;
c = 0.02; % Hệ số giảm chấn
K_0 = 1; % gain
%--------------------------------------------------------------------------
% Chọn n,r, thời gian mô phỏng
tmax = 15;
n = 9; r = 15000;
%--------------------------------------------------------------------------
dt = tmax/(r - 1);
ds = L/(n - 1); 
F1 = zeros(1,r);
F2 = zeros(1,r);
F2(1:2000) = mh*g;
h = zeros(1,r); h(1:3) = ds; % Vị trí xe nâng
w3 = zeros(1,r); % Độ lắc tại vị trí xe nâng
l = 2; % vị trí ban đầu của xe nâng
w = zeros(n,r);
%--------------------------------------------------------------------------
% Chọn bộ điều khiển + SP 
T_set = 3;
sp1 = 1; % Giá trị đặt xe con
sp2 = zeros(1,r); % Giá trị đặt xe nâng
sp2(1:r/2) = 0.3;
sp2(r/2:r) = 0.3;
% Cờ thay đổi giữa các bộ điều khiển
% flag = 0: Check 
% flag = 1: Robust -PD
% flag = 2: Barrier
flag = 0;
w_change = 0;
l_change = 0 ;
e = zeros(1,r);
eh = zeros(1,r);
%--------------------------------------------------------------------------
% Chọn thông số cho bộ quan sát trạng thái 
sp3 = 2;
b0 = 1/(13.1);
b0_n = 1/(0.86);
T_sam = 0.001;
T_sam_n = 0.001;
s_cl = -6/T_set;
s_eso = 9*s_cl;
z_eso = exp(s_eso*T_sam);
Kd = -2*s_cl;                                
Kp = (s_cl)^2;
% l1 = -3*s_eso;
% l2 = 3*(s_eso)^2;
% l3 = -(s_eso)^3;
l1 = 1 - (z_eso)^3;
l2 = 3/(2*T_sam)*(1 - z_eso)^2*(1 + z_eso);
l3 = 1/T_sam^2*(1 - z_eso)^3;
% Các biến trạng thái
x1 = zeros(1,r);
x2 = zeros(1,r);
x3 = zeros(1,r);
x1_n = zeros(1,r);
x2_n = zeros(1,r);
x3_n = zeros(1,r);
%--------------------------------------------------------------------------
% Chọn bộ thông số IS
% ETM4
sp1IS = zeros(1,r);
f = 6.13415;
% f=1.86692;
w_0 = 2*3.14*6.13415;
k_is = exp(c*3.14/sqrt(1 - c^2));
m_is = 4;
i_is = (1 + m_is)*k_is^2/(k_is^2 + (1 + m_is)*(k_is^(4/3) + k_is^(2/3)) + m_is);
a1 = i_is/(1 + m_is);
a2 = i_is/k_is^(2/3);
a3 = i_is/k_is^(4/3);
a4 = m_is*i_is/((1 + m_is)*k_is^2);
t1 = 0;
t2 = (2*3.14/3)/(2*3.14*f);
t3 = (4*3.14/3)/(2*3.14*f);
t4 = (2*3.14)/(2*3.14*f);
%--------------------------------------------------------------------------
%                          MÔ PHỎNG
%--------------------------------------------------------------------------
for j = 3:(r - 1)
    wsss0 = (w(3,j) - 2*w(2,j) + w(1,j))/(2*ds^3);
    w(1,j + 1) = (F1(j + 1) - EI*wsss0)*(dt^2/mw) + 2*w(1,j) - w(1,j - 1);
    % w(1,j + 1) = (F1(j + 1) - EI*wsss0+2*sin(j*pi/4))*(dt^2/mw) + 2*w(1,j) - w(1,j - 1);
    
    for i = 3:(n - 2)
        S2 = (-EI*dt^2)/(ds^4*pA);
        S21 = (-EI*dt^2)/ds^3*mh;
        C = -c*dt^2*(w(i,j) - w(i,j - 1))/dt; % Lực ma sát
        wssss = w(i + 2,j) - 4*w(i + 1,j) + 6*w(i,j) - 4*w(i - 1,j) + w(i - 2,j);
        ddtx3 = (w3(j) - 2*w3(j - 1) + w3(j - 2))/(dt^2);
        dyx2 = (w(l + 1,j - 1) - w(l - 1,j - 1))/2*ds;
        h(j + 1) = (F2(j + 1) - mh*g - mh*ddtx3*dyx2)*dt^2/mh + 2*h(j) - h(j - 1);

        % Cập nhật vị trí xe nâng
        if h(j + 1) < ds
            l = 2;
            h(j + 1) = ds;
        end
        if h(j + 1) > L - ds
            l = n - 1;
            h(j + 1) = L - ds;
        end
        if h(j + 1) > ds && h(j + 1) <= L - ds
            l = ceil(h(j + 1)/ds);
        end

        % Độ lắc của thanh
        if i == l && l > 2 
            w(i,j + 1) = S21*wssss + C + 2*w(i,j) - w(i,j - 1);
            w3(j + 1) = w(i,j + 1) - w(1,j + 1);
        end
        if i ~= l && l >2 
            w(i,j + 1) = S2*wssss + C + 2*w(i,j) - w(i,j - 1);
            w3(j + 1) = w(i,j + 1) - w(1,j + 1);
        end
        if l == 2
            w(i,j + 1) = S2*wssss + C + 2*w(i,j) - w(i,j - 1);
            w3(j + 1) = 0;
        end
        S3 = (EI*(dt^2))/(mk*2*ds^3);
        w(2,j + 1) = w(1,j + 1);
        wsssl = (-2*w(n,j) + 3*w(n - 1,j) - w(n - 2,j));
        w(n,j + 1) = 2*w(n,j) - w(n,j - 1) + S3*wsssl;
        w(n - 1,j + 1) = (w(n,j + 1) + w(n - 2,j + 1))/2;
%--------------------------------------------------------------------------
%                               Bộ điều khiển
%--------------------------------------------------------------------------
        if flag == 0
            F1(1:r/2) = 10;
            F1(r/2:r) = 0;
            F2(1:r) = 20;
            check2 = zeros(1,r);
        end
%--------------------------------------------------------------------------
%                                PD + Robust
%--------------------------------------------------------------------------
        if flag == 1 
            % Thông số robust - PD
            k1 = -5; kx = 20;
            kn = 1;
%--------------------------------------------------------------------------            
            e(j + 1) = (w(1,j + 1) - w(1,j));
            de = (e(j + 1) - e(j))/dt;
            F1(j + 2) = k1*(w(1,j + 1) - sp1) - kx*de;
            % if F1(j + 2) > 50
            %     F1(j + 2) = 50;
            % end
            % if F1(j + 2) < -50
            %     F1(j + 2) = -50;
            % end
            S2 = (w(l + 1,j + 1) - 2*w(l,j + 1) + w(l - 1,j + 1))/(ds^3);
            dx2 = (w(l + 1,j + 1) - w(l - 1,j + 1))/2*ds;
            dtx2 = (h(j + 1) - h(j))/dt;
            F2(j + 2) = mh*g + EI*S2*dx2 - (h(j + 1) - sp2(j)) - kn*dtx2;
        end
%--------------------------------------------------------------------------
%                                 Barrier
%--------------------------------------------------------------------------
        if flag == 2    
            % Thông số 
            k1 = 10; k2 = 0; k3 = 20; k0 = 1;
            kc = 1;
            z = w(l,j + 1) - w(1,j + 1);
            kmax = 0.5;
%--------------------------------------------------------------------------
            dwl = (w(1,j + 1) - w(1,j))/dt;
            S1 = kc/((kmax^2 - z^2)*2.3) + k2;
            F1(j + 2) = -k1*(w(1,j + 1) - sp1) - k3*dwl - S1*kmax*k0*dwl - S1*kmax;
            
            kn = 1;
            S2 = (w(l + 1,j + 1) - 2*w(l,j + 1) + w(l - 1,j + 1))/(ds^3);
            dx2 = (w(l + 1,j + 1) - w(l - 1,j + 1))/2*ds;
            dtx2 = (h(j + 1) - h(j))/dt;
            F2(j + 2) = mh*g + EI*S2*dx2 - (h(j + 1) - sp2(j)) - kn*dtx2; 
        end
    end
end
%--------------------------------------------------------------------------
%                               Vẽ đồ thị
%--------------------------------------------------------------------------
x = 0:dt:tmax;
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

figure(3)
%Phân tích phổ
for i=1:n
Fs=1/dt;
time=0:dt:(r-2)*dt-dt;
l1=r;
fft_w=fft(w(i,:)-w(1,:),l1)*(2/l1);
abs_w=abs(fft_w);
freq=0:(1/time(end)):Fs/2-(1/time(end));
subplot(ceil(n/2), 2, i);
plot(freq,abs_w(1:length(freq)),'LineWidth',0.8);
hold on
grid on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 9);
xlabel('Frequency(Hz)');
ylabel('Amplitude(m)');
% title('Pho dao dong');
title(['FFT tại vị trí x = ', num2str(ds*(i-1)), ' m']);
end

