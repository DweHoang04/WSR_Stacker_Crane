clc;
close all;
clear;
%--------------------------------------------------------------------------
%    Mô phỏng chuyển động của dầm Euler - Bernoulli trong mô hình SMC
%--------------------------------------------------------------------------

% Thiết lập các thông số cho dầm
L = 0.63; EI = 0.754; rho_A = .297;
% Các vật nặng 
mw = 13.1; mk = 0.04; mh = 0.86; g = 9.81;

% Thiết lập thông số không gian và thời gian
n = 5; r = 5000;
tmax = 15;
delta_Y = L/(n - 1); % Bước không gian
delta_t = tmax/(r - 1); % Bước thời gian

% Khởi tạo ma trận để lưu giá trị
w = zeros(n,r);
x3 = zeros(1,r); % Độ lắc của thanh tại x2 (x3)
x2 = 2; 
dx2dt_2 = zeros(1,r);
dx2dt_2(1:2) = delta_Y;

% Lực tác động vào xe con
F1 = zeros(1,r);
F1(1:r) = 10;

% Lực tác động vào xe nâng
F2 = zeros(1,r);
F2(1:1000) = 10;

%--------------------------------------------------------------------------
for j = 3:(r - 1)
    wyyy0 = (w(3,j) - 2*w(2,j) + w(1,j))/(2*delta_Y^3);
    S2 = (F1(j + 1) - EI*wyyy0)/mw;
    w(1,j + 1) = 2*w(1,j) - w(1,j - 1) + delta_t^2*S2; % 5b
    for i = 3:(n - 2)
        % Đạo hàm theo Y
        wyyyy = (w(i + 2,j) - 4*w(i + 1,j) + 6*w(i,j) - 4*w(i - 1,j) + w(i - 2,j));
        S1 = (-EI/(rho_A*delta_Y^4))*wyyyy; 
        S4 = (-EI/(mh*delta_Y^3))*wyyyy;

        % Chuyển động của xe nâng
        dx3dt_2 = (x3(j + 1) - 2*x3(j) + x3(j - 1))/(2*delta_t^2);
        wyx2 = (w(x2 + 1,j - 1) - w(x2 - 1,j - 1))/(2*delta_Y^2);
        dx2dt_2(j + 1) = 2*dx2dt_2(j) - dx2dt_2(j - 1) + (F2(j + 1) - mh*g - mh*dx3dt_2*wyx2)/mh*delta_t^2; % 5c

        % Cập nhật vị trí x2
        if dx2dt_2(j + 1) < delta_Y
            x2 = 2;
            dx2dt_2(j + 1) = delta_Y;
        elseif dx2dt_2(j + 1) > L - delta_Y
            x2 = n - 1;
            dx2dt_2(j + 1) = L - delta_Y;
        else
            x2 = ceil(dx2dt_2(j + 1)/delta_Y);
        end

        % Độ lắc của thanh
        if dx2dt_2(j + 1) ~= x2
            w(i,j + 1) = 2*w(i,j) - w(i,j - 1) + delta_t^2*S1; % 5a
            x3(j + 1) = w(i,j + 1) - w(1,j + 1);
        else
            w(i,j + 1) = 2*w(i,j) - w(i,j - 1) + delta_t^2*S4; % 5d
            x3(j + 1) = w(i,j + 1) - w(1,j + 1);
        end
    end
    wyyyl = (-2*w(n,j) + 3*w(n - 1,j) - w(n - 2,j))/(2*delta_Y^3);
    S3 = (EI/mk)*wyyyl;
    w(2,j + 1) = w(1,j + 1);
    w(n,j + 1) = 2*w(n,j) - w(n,j - 1) + delta_t^2*S3; % 5e
    w(n - 1,j + 1) = (w(n,j + 1) + w(n - 2,j + 1))/2;
end
%--------------------------------------------------------------------------

% Khởi tạo giá trị để mô phỏng
t_tr = linspace(0,tmax,r);

% Vẽ đồ thị
subplot(2,2,1);
grid on;
hold on;
plot(t_tr,w(1,:),'b','LineWidth',1.5);
title({'Vị trí của xe con'});
ylabel('Vị trí xe con (m)','FontSize',12);
xlabel('Thời gian (s)','FontSize',12);

subplot(2,2,2);
grid on;
hold on;
plot(t_tr,w(n,:) - w(1,:),'g','LineWidth',1.5); % Vị trí tương đối của mk so với xe con
title({'Độ lắc đỉnh thanh'});
ylabel('Độ lắc so với vị trí cân bằng (m)','FontSize',12);
xlabel('Thời gian (s)','FontSize',12);

subplot(2,2,3);
grid on;
hold on;
plot(t_tr,dx2dt_2,'r','LineWidth',1.5);
title({'Vị trí xe nâng'});
ylabel('Vị trí xe nâng (m)','FontSize',12);
xlabel('Thời gian (s)','FontSize',12);

subplot(2,2,4);
grid on;
hold on;
plot(t_tr,x3,'Color',[1 0.5 0],'LineWidth',1.5);
title({'Độ lắc xe nâng'});
ylabel('Độ lắc xe nâng (m)','FontSize',12);
xlabel('Thời gian (s)','FontSize',12);