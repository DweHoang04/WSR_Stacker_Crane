clc;
close all;
clear;
%--------------------------------------------------------------------------
%                         Cài đặt các thông số
%--------------------------------------------------------------------------

% Thiết lập các thông số cho dầm
L = 0.63; EI = 0.754; rho_A = 0.297;
% Các vật nặng 
mw = 13.1; mk = 0.04; mh = 0.86; g = 9.81;

% Thiết lập thông số không gian và thời gian
n = 9; r = 20000;
tmax = 20;
delta_Y = L/(n - 1); % Bước không gian
delta_t = tmax/(r - 1); % Bước thời gian

% Khởi tạo ma trận để lưu giá trị
w = zeros(n,r);
x3 = zeros(1,r); % Độ lắc của thanh tại x2 (x3)
x2 = 2; 
dx2dt_2 = zeros(1,r);
dx2dt_2(1:3) = delta_Y;

% Lực tác động vào xe con
F1 = zeros(1,r);
F1(1:r) = 10;

% Lực tác động vào xe nâng
F2 = zeros(1,r);
F2(1:r) = 3;

%--------------------------------------------------------------------------
%                        Thông số bộ điều khiển ADRC
%--------------------------------------------------------------------------
T_set = 2; % Thời gian xác lập (s)
T_sample = 0.001; % Chu kỳ trích mẫu
s_CL = -6/T_set; % Điểm cực hàm truyền hệ kín
% Hệ số tỉ lệ
Kp = s_CL^2; Kd = -2*s_CL; K_ESO = 10; b01 = 1/mw; b02 = 1/mh;

% Bộ quan sát trạng thái ESO
s_ESO = K_ESO*s_CL; % Điểm cực 
z_ESO = exp(s_ESO*T_sample);

l1 = 1 - z_ESO^3; 
l2 = (3*(1 - z_ESO)^2*(1 + z_ESO))/(2*T_sample);
l3 = (1 - z_ESO)^3/(T_sample^2);

Ad = [1 T_sample (T_sample^2)/2; 0 1 T_sample; 0 0 1];
Bd = [b01*(T_sample^2)/2; b01*T_sample; 0];
Cd = [1 0 0];
Lc = [l1; l2; l3];

A_ESO = Ad - Lc*Cd*Ad;
B_ESO = Bd - Lc*Cd*Bd;

% Biến trạng thái xk = [x1k; x2k; x3k]
xd = zeros(3,r); 
xl = zeros(3,r);

x1_set = 1; % Giá trị đặt xe con
x2_set = 0.4; % Giá trị đặt xe nâng

%--------------------------------------------------------------------------
%    Mô phỏng chuyển động của dầm Euler - Bernoulli trong mô hình SMC
%--------------------------------------------------------------------------
for j = 3:(r - 1)
    wyyy0 = (w(3,j) - 2*w(2,j) + w(1,j))/(2*delta_Y^3);
    S2 = (F1(j + 1) - EI*wyyy0)/mw;
    w(1,j + 1) = 2*w(1,j) - w(1,j - 1) + delta_t^2*S2; % 5b
    for i = 3:(n - 2)
        % Đạo hàm theo Y
        wyyyy = w(i + 2,j) - 4*w(i + 1,j) + 6*w(i,j) - 4*w(i - 1,j) + w(i - 2,j);
        S1 = (-EI/(rho_A*delta_Y^4))*wyyyy;
        S4 = (-EI/(mh*delta_Y^3))*wyyyy;
        w(i,j + 1) = 2*w(i,j) - w(i,j - 1) + delta_t^2*S1; % 5a

        % Chuyển động của xe nâng
        dx3dt_2 = (x3(j) - 2*x3(j - 1) + x3(j - 2))/(2*delta_t^2);
        wyx2 = (w(x2 + 1,j - 1) - w(x2 - 1,j - 1))/(2*delta_Y^2);
        dx2dt_2(j + 1) = 2*dx2dt_2(j) - dx2dt_2(j - 1) + ((F2(j + 1) - mh*g - mh*dx3dt_2*wyx2)*delta_t^2)/mh; % 5c

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
%--------------------------------------------------------------------------
%                 Bộ điều khiển ADRC cho xe con và xe nâng
%--------------------------------------------------------------------------
    % Xe con    
    xd(:,j + 2) = A_ESO*xd(:,j + 1) + B_ESO*F1(j + 1) + Lc*w(1,j + 1);
    F1(j + 2) = (Kp*(x1_set - xd(1,j + 2)) - Kd*xd(2,j + 2) - xd(3,j + 2))/b01;

    % Xe nâng
    xl(:,j + 2) = A_ESO*xl(:,j + 1) + B_ESO*(F2(j + 1) - mh*g) + Lc*dx2dt_2(j + 1);
    F2(j + 2) = (Kp*(x2_set - xl(1,j + 2)) - Kd*xl(2,j + 2) - xl(3,j + 2))/b02 + mh*g;
end
%--------------------------------------------------------------------------

% Khởi tạo giá trị để mô phỏng
t_tr = 0:delta_t:tmax;

% Vẽ đồ thị
figure(1);
subplot(2,2,1);
grid on;
hold on;
plot(t_tr,w(1,:),'g','LineWidth',1.5);
title({'Vị trí của xe con'});
ylabel('Vị trí xe con (m)','FontSize',12);
xlabel('Thời gian (s)','FontSize',12);

subplot(2,2,2);
grid on;
hold on;
plot(t_tr,w(n,:) - w(1,:),'b','LineWidth',1.5); % Vị trí tương đối của mk so với xe con
title({'Độ lắc đỉnh thanh'});
ylabel('Biên độ dao động (m)','FontSize',12);
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

figure(2);
subplot(1,2,1);
grid on;
hold on;
plot(t_tr,F1(1:r),'LineWidth',1.5);
ylabel('Độ lớn lực F1 (N)','FontSize',12);
xlabel('Thời gian (s)','FontSize',12);

subplot(1,2,2);
grid on;
hold on;
plot(t_tr,F2(1:r),'LineWidth',1.5);
ylabel('Độ lớn lực F2 (N)','FontSize',12);
xlabel('Thời gian (s)','FontSize',12);
