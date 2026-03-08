clc;
close all;
clear;
%--------------------------------------------------------------------------
%    Mô phỏng chuyển động của dầm Euler - Bernoulli trong mô hình SMC
%--------------------------------------------------------------------------

% Thiết lập các thông số cho dầm
L = 1; EI = 14.97; rho_A = 2.1; % EI = 9.6 thì dao động lên được 5*10^-3
% Các vật nặng 
mw = 15; mk = 0.2; mh = 0.9; g = 9.81;

% Thiết lập thông số không gian và thời gian
n = 9; r = 10000;
tmax = 10;
delta_Y = L/(n - 1); % Bước không gian
delta_t = tmax/(r - 1); % Bước thời gian

% Khởi tạo ma trận để lưu giá trị
w = zeros(n,r);         % Độ lắc cả thanh
x3 = zeros(1,r);        % Độ lắc của thanh tại x2 (x3)
x2 = 2;                 % Vị trí ban đầu của xe nâng
wx2 = zeros(1,r);       % Độ lắc tại vị trí xe nâng
wx2(1:3) = delta_Y;

% Lực tác động vào xe con
F1 = zeros(1,r);
F1(1:1000) = 10;
F1(1001:2000) = -10;

% % Lực tác động vào xe nâng
F2 = zeros(1,r);
F2(1:2000) = 9;

% F2(1:r)=mh*g;
% F1(1:1500)=10;
% F2(1:1500)=mh*g+10;
% F1(1501:3000)=-10;

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
        wx2(j + 1) = 2*wx2(j) - wx2(j - 1) + (F2(j + 1) - mh*g - mh*dx3dt_2*wyx2)/mh*delta_t^2; % 5c

        % Cập nhật vị trí x2
        if wx2(j + 1) < delta_Y
            x2 = 2;
            wx2(j + 1) = delta_Y;
        elseif wx2(j + 1) > L - delta_Y
            x2 = n - 1;
            wx2(j + 1) = L - delta_Y;
        else
            x2 = ceil(wx2(j + 1)/delta_Y);
        end

        % Độ lắc của thanh
        if wx2(j + 1) ~= x2*delta_Y
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
t_tr = 0:delta_t:tmax;

% Vẽ đồ thị
figure(1);
grid on; hold on;
plot(t_tr, w(1,:), 'Color',[0.07,0.62,1.00], 'LineWidth',4);
title('\textbf{Driving Unit Position}', 'Interpreter','latex', 'FontSize', 30);
ylabel('$x_1(t)\ \mathrm{(m)}$', 'Interpreter','latex', 'FontSize', 20);
xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 20);
axis([0 10 0 1]);
% Plot points
plot(1,  0.3483,'o','MarkerSize',7,'MarkerEdgeColor','black','LineWidth',2);
plot(2,  0.6960,'o','MarkerSize',7,'MarkerEdgeColor','black','LineWidth',2);
plot(10, 0.6792,'o','MarkerSize',7,'MarkerEdgeColor','black','LineWidth',2);
% Coordinate labels
text(1+0.1,  0.3483, '$(1,\;0.3483)$', ...
    'Interpreter','latex', 'FontSize',20);
text(2-0.43,  0.6960+0.03, '$(2,\;0.6960)$', ...
    'Interpreter','latex', 'FontSize',20);
text(10-1.2, 0.6792+0.03, '$(10,\;0.6792)$', ...
    'Interpreter','latex', 'FontSize',20);
%-----------------------------------------------------------------------
figure(2)
grid on; hold on;
plot(t_tr, wx2, 'Color',[1 0.5 0], 'LineWidth',4);
title('\textbf{Lifting Unit Position}', 'Interpreter','latex', 'FontSize', 30);
ylabel('$x_2(t)\ \mathrm{(m)}$', 'Interpreter','latex', 'FontSize', 20);
xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 20);
axis([0 10 0 1]);
plot([0 2.431], [0.125 0.125], 'k--', 'LineWidth',3);
plot(1.462, 0.875,'o','MarkerSize',7,'MarkerEdgeColor','black','LineWidth',2);
plot(2.431, 0.125,'o','MarkerSize',7,'MarkerEdgeColor','black','LineWidth',2);
text(1.462-0.58, 0.875+0.04, '$(1.462,\;0.875)$', ...
    'Interpreter','latex', 'FontSize',20);
text(2.431-0.58, 0.125-0.04, '$(2.431,\;0.125)$', ...
    'Interpreter','latex', 'FontSize',20);

figure(3)
grid on; hold on;
y_beam = w(n,:) - w(1,:);
plot(t_tr, y_beam, 'Color',[0.07,0.62,1.00], 'LineWidth',2);
title('\textbf{Beam''s Top Vibration}', 'Interpreter','latex', 'FontSize', 30);
ylabel('$x_4(t)\ \mathrm{(m)}$', 'Interpreter','latex', 'FontSize', 20);
xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 20);
axis([0 10 -3.5e-3 3.5e-3]);

% Zoom inset
zoom_x_start = 3;
zoom_x_end   = 4;
zoom_y_start = -3.5e-3;
zoom_y_end   =  3.5e-3;
rectangle('Position', ...
    [zoom_x_start, zoom_y_start, ...
     zoom_x_end-zoom_x_start, zoom_y_end-zoom_y_start], ...
    'EdgeColor','red','LineWidth',2);
axes('Position',[0.60 0.47 0.35 0.4]);
box on; hold on; grid on;
plot(t_tr, y_beam, 'Color',[0.07,0.62,1.00], 'LineWidth',1.2);
xlim([zoom_x_start zoom_x_end]);
ylim([zoom_y_start zoom_y_end]);
set(gca,'FontSize',8);
annotation('arrow', [0.44 0.6], [0.5 0.67], ...
    'Color','red','LineWidth',2,'HeadStyle','vback2');

y_max =  2.9e-3;
y_min = -2.9e-3;
plot([3 4], [y_max y_max], 'r--', 'LineWidth',2);
plot([3 4], [y_min y_min], 'r--', 'LineWidth',2);
x_arrow = 3.5;
line([x_arrow x_arrow], [y_min y_max], 'Color','r','LineWidth',1.5);
plot(x_arrow, 2.8e-3, '^r', 'MarkerFaceColor','r');
plot(x_arrow, -2.8e-3, 'vr', 'MarkerFaceColor','r');
text(x_arrow-0.125, 0, '$5.6\times10^{-3}\ \mathrm{m}$', ...
    'Interpreter','latex','FontSize',20, ...
    'BackgroundColor','white','EdgeColor','red');
title('\textbf{Zoom View}','Interpreter','latex','FontSize',20);

%-----------------------------------------------------------------------
% Dao động cả thanh
figure('Position', [100 100 600 400]);
xc = w - w(1,:);
xc1 = xc';
[X, Y] = meshgrid(0:delta_t:tmax,0:delta_Y:L);
meshc(Y, X, xc);
title('\textbf{3D View Beam Vibration}', 'Interpreter','latex', 'FontSize', 30);
xlabel('$Y(t)\ \mathrm{(m)}$', 'Interpreter','latex', 'FontSize', 20);
ylabel('Time (s)', 'Interpreter','latex', 'FontSize', 20);
zlabel('$w(t)\ \mathrm{(m)}$', 'Interpreter','latex', 'FontSize', 20);
view(65,10);
axis([0 1 0 10 -3.5*10^-3 3.5*10^-3]);
yticks([0 2 4 6 8 10]);