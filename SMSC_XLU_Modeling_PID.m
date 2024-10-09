clc;
close all;
clear;
%--------------------------------------------------------------------------
%    Mô phỏng chuyển động của dầm Euler - Bernoulli trong mô hình SMC
%--------------------------------------------------------------------------

% Thiết lập các thông số cho dầm
L = 0.63; EI = 0.754; rho_A = 0.297;
% Các vật nặng 
mw = 13.1; mk = 0.04; mh = 0.86; g = 9.81;

% Thiết lập thông số không gian và thời gian
n = 9; r = 20000;
tmax = 30;
delta_Y = L/(n - 1); % Bước không gian
delta_t = tmax/(r - 1); % Bước thời gian

% Bộ thông số PID
kp1 = 0.1; ki1 = 3; kd1 = 9; % Hệ số 3 khâu PID cho xe con
I = 0; % Khởi tạo giá trị tính tích phân
sp1 = 1; % Giá trị đặt cho xe con
e1 = zeros(1,r); % Sai số vị trí xe con

% Khởi tạo ma trận để lưu giá trị
w = zeros(n,r);

% Lực tác động vào xe con
F1 = zeros(1,r);
F1(1:r/2) = 1;

%--------------------------------------------------------------------------
for j = 3:(r - 1)
    wyyy0 = (w(3,j) - 2*w(2,j) + w(1,j))/(2*delta_Y^3);
    S2 = (F1(j + 1) - EI*wyyy0)/mw;
    w(1,j + 1) = 2*w(1,j) - w(1,j - 1) + delta_t^2*S2; % 5b
    for i = 3:(n - 2)
        % Đạo hàm theo Y
        wyyyy = (w(i + 2,j) - 4*w(i + 1,j) + 6*w(i,j) - 4*w(i - 1,j) + w(i - 2,j));
        S1 = (-EI/(rho_A*delta_Y^4))*wyyyy;
        w(i,j + 1) = 2*w(i,j) - w(i,j - 1) + delta_t^2*S1; % 5a
    end
    wyyyl = (-2*w(n,j) + 3*w(n - 1,j) - w(n - 2,j))/(2*delta_Y^3);
    S3 = (EI/mk)*wyyyl;
    w(2,j + 1) = w(1,j + 1);
    w(n,j + 1) = 2*w(n,j) - w(n,j - 1) + delta_t^2*S3; % 5e
    w(n - 1,j + 1) = (w(n,j + 1) + w(n - 2,j + 1))/2;
%--------------------------------------------------------------------------
%                       Bộ điều khiển PD cho xe con
%--------------------------------------------------------------------------
    e1(j + 1) = sp1 - w(1,j);
    I = I + e1(j + 1)*delta_t;
    D = (e1(j + 1) - e1(j))/delta_t;
    F1(j + 2) = kp1*e1(j + 1) + kd1*D;
end
%--------------------------------------------------------------------------

% Khởi tạo giá trị để mô phỏng
t_tr = linspace(0,tmax,r);

% Vẽ đồ thị
figure(1)
subplot(2,1,1);
grid on;
hold on;
plot(t_tr,w(1,:),'b','LineWidth',1.5);
title({'Vị trí của xe con'});
ylabel('Vị trí xe con (m)','FontSize',12);
xlabel('Thời gian (s)','FontSize',12);

subplot(2,1,2);
grid on;
hold on;
plot(t_tr,w(n,:) - w(1,:),'g','LineWidth',1.5); % Vị trí tương đối của mk so với xe con
title({'Độ lắc đỉnh thanh'});
ylabel('Độ lắc so với vị trí cân bằng (m)','FontSize',12);
xlabel('Thời gian (s)','FontSize',12);

figure(2)
grid on;
hold on;
plot(t_tr,F1(1:r),'g','LineWidth',1.5);
title({'Lực F1'});
ylabel('Lực F1','FontSize',12);
xlabel('Thời gian (s)','FontSize',12);