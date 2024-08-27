clc;
close all;
clear;
%--------------------------------------------------------------------------
%    Mô phỏng chuyển động của dầm Euler - Bernoulli trong mô hình SMC
%--------------------------------------------------------------------------

% Thiết lập các thông số cho dầm
L = 1.5; EI = 14.97; rho_A = .21;
% Các vật nặng 
mw = 13.1; mk = 0.19;

% Thiết lập thông số không gian và thời gian
n = 3; r = 3000;
tmax = 15;
delta_Y = L/(n - 1); % Bước không gian
delta_t = tmax/(r - 1); % Bước thời gian

% Khởi tạo ma trận để lưu giá trị
w = zeros(n,r);

w_3D_free = zeros(r,n);

% Lực tác động vào xe con
F1 = zeros(1,r);
F1(1:r) = 10;

for j = 2:(r - 1)
    for i = 3:(n - 2)
        % Đạo hàm theo Y
        wyyyy = (w(i + 2,j) - 4*w(i + 1,j) + 6*w(i,j) - 4*w(i - 1,j) + w(i - 2,j))/delta_Y^4;

        S1 = (-EI/rho_A)*wyyyy;
        
        % Chuyển động của dầm Euler - Bernoulli
        w(i,j + 1) = 2*w(i,j) - w(i,j - 1) + delta_t^2*S1;
    end

    % Đạo hàm theo Y tại chân và đỉnh của dầm
    wyyy0 = (w(3,j) - 2*w(2,j) + w(1,j))/(2*delta_Y^3);
    wyyyl = (-2*w(n,j) + 3*w(n - 1,j) - w(n - 2,j))/(2*delta_Y^3);

    S2 = (F1(j) - EI*wyyy0)/mw;
    S3 = (EI/mk)*wyyyl;

    % Dao động tại chân, chính giữa và đỉnh của dầm
    w(1,j + 1) = 2*w(1,j) - w(1,j - 1) + delta_t^2*S2;
    %w(2,j + 1) = w(1,j + 1);
    w(ceil(n/2),j + 1) = (w(1,j + 1) + w(n,j + 1))/2;
    w(n,j + 1) = 2*w(n,j) - w(n,j - 1) + delta_t^2*S3;
    
    w(n - 1,j + 1) =  (w(n,j + 1) + w(n - 2,j + 1))/2;

    w_3D_free(1 + j,:) = w(:,j)';
end

w_3D_free(1,:) = w(:,1)';

% Khởi tạo giá trị để mô phỏng
w0 = linspace(0,1,n);
t_tr = linspace(0,tmax,r);

% Vẽ biểu đồ
figure(1);
grid on;
hold on;
surf(w0,t_tr,w_3D_free); view(45,30);
title({'Dao động của thanh khi không có điều khiển'});
ylabel('t','FontSize',12);
xlabel('y','FontSize',12);
zlabel('w(Y,t)','FontSize',12);

figure(2)
grid on;
hold on;
plot(t_tr,w(1,:));
title({'Vị trí của xe con'});
ylabel('x1','FontSize',12);
xlabel('t','FontSize',12);

figure(3)
grid on;
hold on;
plot(t_tr,w(n,:) - w(1,:)); % Vị trí tương đối của thanh so với xe con
title({'Vị trí của đỉnh thanh'});
ylabel('x4','FontSize',12);
xlabel('t','FontSize',12);