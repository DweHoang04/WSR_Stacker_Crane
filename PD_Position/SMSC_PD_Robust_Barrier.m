clc;
close all;
clear;

% Cờ để thay đổi giữa các bộ điều khiển
% flag = 0: Không điều khiển; flag = 1: PD + Robust; flag = 2: Barrier
flag = 2;

%--------------------------------------------------------------------------
%    Mô phỏng chuyển động của dầm Euler - Bernoulli trong mô hình SMC
%--------------------------------------------------------------------------

% Thiết lập các thông số cho dầm
L = 0.63; EI = 0.754; rho_A = 0.297;
% Các vật nặng 
mw = 13.1; mk = 0.04; mh = 0.86; g = 9.81;

% Thiết lập thông số không gian và thời gian
n = 10; r = 20000;
tmax = 15;
delta_Y = L/(n - 1); % Bước không gian
delta_t = tmax/(r - 1); % Bước thời gian

%--------------------------------------------------------------------------
%                   Thông số bộ điều khiển PD + Robust
%--------------------------------------------------------------------------
k1 = 10; kx = 20; kn = 1.5;
x1_set = 1; % Giá trị đặt xe con
x2_set = 0.3; % Giá trị đặt xe nâng
%--------------------------------------------------------------------------
%                     Thông số bộ điều khiển Barrier
%--------------------------------------------------------------------------
k1_barrier = 20; kmax = 0.004; k2 = 0.01; k3 = 30; 
k0 = 10; kc = 0.001; kn_barrier = 1.5;
%--------------------------------------------------------------------------

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
F2(1:r) = 0;

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
        w(i,j + 1) = 2*w(i,j) - w(i,j - 1) + delta_t^2*S1; % 5a

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
%--------------------------------------------------------------------------
%                        Bộ điều khiển PD + Robust
%--------------------------------------------------------------------------
    if flag == 1
        % e(j + 1) = (w(1,j + 1) - w(1,j));
        % de = (e(j + 1) - e(j))/delta_Y;
        dx1 = (w(1,j + 1) - w(1,j))/delta_t;
        F1(j + 2) = -k1*(w(1,j + 1) - x1_set) - kx*dx1;
        % if F1(j + 2)>50
        %     F1(j + 2) = 50;
        % end
        % if F1(j + 2) < -5
        %     F1(j + 2) = -50;
        % end
        S2 = (w(x2 + 1,j + 1) - 2*w(x2,j + 1) + w(x2 - 1,j + 1))/(delta_Y^3);
        dx2 = (w(x2 + 1,j + 1) - w(x2 - 1,j + 1))/2*delta_Y;
        dtx2 = (dx2dt_2(j + 1) - dx2dt_2(j))/delta_t;
        F2(j + 2) = mh*g + EI*S2*dx2 - (dx2dt_2(j + 1) - x2_set) - kn*dtx2;
    end
%--------------------------------------------------------------------------
%                       Bộ điều khiển PD + Barrier
%--------------------------------------------------------------------------
    if flag == 2
        z = w(x2,j + 1);
        dwl = (w(1,j + 1) - w(1,j))/delta_t;
        S1 = kc/((kmax^2 - z^2)*2.3) + k2;
        F1(j + 2) = -k1_barrier*(w(1,j + 1) - x1_set) - k3*dwl ...
                    - S1*kmax*k0*dwl - S1*kmax;
        
        S2 = (w(x2 + 1,j + 1) - 2*w(x2,j + 1) + w(x2 - 1,j + 1))/(delta_Y^3);
        dx2 = (w(x2 + 1,j + 1) - w(x2 - 1,j + 1))/2*delta_Y;
        dtx2 = (dx2dt_2(j + 1) - dx2dt_2(j))/delta_t;
        F2(j + 2) = mh*g + EI*S2*dx2 - (dx2dt_2(j+1) - x2_set) - kn*dtx2;
    end
end
%--------------------------------------------------------------------------

% Khởi tạo giá trị để mô phỏng
t_tr = linspace(0,tmax,r);

% Vẽ đồ thị
figure(1) % Đồ thị vị trí và dao động
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

% Đồ thị dao động 3D
% n1 = 2; r1 = 1000;
% wg = zeros(n1,r1);
% rn1 = n/n1; rr1 = r/r1;
% for i = 1:n1
%     for j = 1:r1
%         wg(i,j) = w(rn1*i - (rn1 - 1), rr1*j - (rr1 - 1));
%     end
% end
% xc = (wg - wg(1,:))*100000;
% xc1 = xc';
% 
% figure('Position', [700 100 600 300]);
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 9);
% [X, Y] = meshgrid(0:(tmax/(r1 - 1)):tmax,0:(L/(n1 - 1)):L);
% meshc(Y,X,xc);
% % ylabel('t(s)');
% % xlabel('Y(m)');
% % zlabel('ω(Y,t)(mm)');
% 
% ylabel('$t$(s)', 'Interpreter', 'latex', 'FontSize',9);
% xlabel('$Y$(m)', 'Interpreter', 'latex', 'FontSize',9);
% zlabel('$\omega(Y,t)$(mm)', 'Interpreter', 'latex', 'FontSize',9);
% view(65,10);
% axis([0 0.6 0 10 -40 40])

