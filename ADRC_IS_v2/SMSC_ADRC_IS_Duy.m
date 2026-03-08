clc;
close all;
clear;

% Cờ để thay đổi giữa ADRC và ADRC + IS
% flag = 1: ADRC; flag = 2: ADRC + IS
flag = 1; 
% Cờ để thay đổi giữa các bộ IS
% flag_is = 1: ZV; flag_is = 2: ZVD; flag = 3: ETM4; flag = 4: TVZV
flag_is = 1;

%--------------------------------------------------------------------------
%                         Cài đặt các thông số
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
% F2(1:r) = mh*g - 0.1;

%--------------------------------------------------------------------------
%                        Thông số bộ điều khiển ADRC
%--------------------------------------------------------------------------
T_set = 3; % Thời gian xác lập (s)
T_sample = delta_t; % Chu kỳ trích mẫu
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
Bd1 = [b01*(T_sample^2)/2; b01*T_sample; 0];
Bd2 = [b02*(T_sample^2)/2; b02*T_sample; 0];

Cd = [1 0 0];
Lc = [l1; l2; l3];

A_ESO = Ad - Lc*Cd*Ad;
B_ESO1 = Bd1 - Lc*Cd*Bd1;
B_ESO2 = Bd2 - Lc*Cd*Bd2;

% Biến trạng thái xk = [x1k; x2k; x3k]
xd = zeros(3,r); 
xl = zeros(3,r);

x1_set = 1; % Giá trị đặt xe con
x2_set = 0.4; % Giá trị đặt xe nâng
%--------------------------------------------------------------------------
%                       Thông số bộ tạo dạng tín hiệu
%--------------------------------------------------------------------------
f = 5.60112; % Tần số dao động riêng của hệ

zeta = 0; K = exp((zeta*pi)/sqrt(1 - zeta^2)); 
mopt = 1; % Giá trị tối ưu cho bộ ETM

% Các thời điểm của các vector xung
t2 = 1/(2*f); t3 = 1/f;

x_IS = zeros(1,r);

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
        dx2dt_2(j + 1) = 2*dx2dt_2(j) - dx2dt_2(j - 1) + ((F2(j + 1) - mh*g ...
                         - mh*dx3dt_2*wyx2)*delta_t^2)/mh; % 5c

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
    if flag == 1
        % Xe con    
        xd(:,j + 2) = A_ESO*xd(:,j + 1) + B_ESO1*F1(j + 1) + Lc*w(1,j + 1);
        F1(j + 2) = (Kp*(x1_set - xd(1,j + 2)) ...
                    - Kd*xd(2,j + 2) - xd(3,j + 2))/b01;
    
        % Xe nâng
        xl(:,j + 2) = A_ESO*xl(:,j + 1) + B_ESO2*(F2(j + 1) - mh*g) ...
                      + Lc*dx2dt_2(j + 1);
        F2(j + 2) = (Kp*(x2_set - xl(1,j + 2)) - Kd*xl(2,j + 2) ...
                    - xl(3,j + 2))/b02 + mh*g;
    end
%-----------------------------------------------------------------------
%                      Bộ tạo dạng tín hiệu IS + ADRC
%-----------------------------------------------------------------------
    if flag == 2
        % IS
        if flag_is == 1 
            % Bộ tạo dạng tín hiệu ZV
            A1 = K/(1 + K); A2 = A1;
            x1_set_is = x1_set*(A1*unit_step(j*delta_t) ...
                        + A2*unit_step(j*delta_t - t2));
        elseif flag_is == 2 
            % Bộ tạo dạng tín hiệu ZVD
            A1 = (K/(1 + K))^2; 
            A2 = (2*K)/(1 + K)^2;
            A3 = 1/(1 + K)^2;
            x1_set_is = x1_set*(A1*unit_step(j*delta_t) + A2*unit_step(j*delta_t - t2) ...
                        + A3*unit_step(j*delta_t - t3));
        elseif flag_is == 3 
            % Bộ tạo dạng tín hiệu ETM4
            t2 = 1/(3*f); t3 = 2/(3*f); t4 = 1/f;
            I = (K^2*(1 + mopt))/(K^2 + (1 + mopt)*(K^(4/3) + K^(2/3)) + mopt);
            A1 = I/(1 + mopt); A2 = I/K^(2/3);
            A3 = I/K^(4/3); A4 = (mopt*I)/(K^2*(1 + mopt));
            x1_set_is = x1_set*(A1*unit_step(j*delta_t) + A2*unit_step(j*delta_t - t2) ...
                        + A3*unit_step(j*delta_t - t3)) + A4*unit_step(j*delta_t - t4);
        end
        x_IS(j) = x1_set_is;
        % ADRC
        % Xe con
        xd(:,j + 2) = A_ESO*xd(:,j + 1) + B_ESO1*F1(j + 1) + Lc*w(1,j + 1);
        F1(j + 2) = (Kp*(x1_set_is - xd(1,j + 2)) - Kd*xd(2,j + 2) - xd(3,j + 2))/b01;
        
        % Xe nâng
        xl(:,j + 2) = A_ESO*xl(:,j + 1) + B_ESO2*(F2(j + 1) - mh*g) + Lc*dx2dt_2(j + 1);
        F2(j + 2) = (Kp*(x2_set - xl(1,j + 2)) - Kd*xl(2,j + 2) - xl(3,j + 2))/b02 + mh*g;
    end
end
%--------------------------------------------------------------------------

% Khởi tạo giá trị để mô phỏng
t_tr = 0:delta_t:tmax;

% Vẽ đồ thị
figure(1);
grid on;
hold on;
% Vị trí xe con
plot(t_tr,w(1,:),'Color',[0.07,0.62,1.00],'LineWidth',4);
title('\textbf{Driving Unit Position}', 'Interpreter','latex', 'FontSize', 30);
ylabel('$x_1(t)\ \mathrm{(m)}$', 'Interpreter','latex', 'FontSize', 20);
xlabel('Time (s)','Interpreter','latex','FontSize', 20);
axis([0 10 0 1.2]);

% Zoom inset
zoom_x_start = 2;
zoom_x_end   = 10;
zoom_y_start = 0.95;
zoom_y_end   =  1.05;
rectangle('Position', ...
    [zoom_x_start, zoom_y_start, ...
     zoom_x_end-zoom_x_start, zoom_y_end-zoom_y_start], ...
    'EdgeColor','red','LineWidth',2);
axes('Position',[0.4 0.2 0.5 0.4]);
box on; hold on; grid on;
plot(t_tr, w(1,:), 'Color',[0.07,0.62,1.00], 'LineWidth',3);
xlim([zoom_x_start zoom_x_end]);
ylim([zoom_y_start zoom_y_end]);
set(gca,'FontSize',8);
p1 = plot(2.915,0.98,'o','MarkerSize',7,'MarkerEdgeColor',[0.5020, 0.0, 0.5020],'LineWidth',1.5,'MarkerIndices',1);
p2 = plot(5.356,1,'o','MarkerSize',7,'MarkerEdgeColor','red','LineWidth',1.5,'MarkerIndices',1);
p3 = plot(5.459,0.9996,'o','MarkerSize',7,'MarkerEdgeColor','black','LineWidth',1.5,'MarkerIndices',1);
p4 = plot(8.929,1.0002,'o','MarkerSize',7,'MarkerEdgeColor','black','LineWidth',1.5,'MarkerIndices',1);

text(2.915+0.3, 0.98, '$(2.915,\;0.98)$', ...
    'Interpreter','latex','FontSize',20, ...
    'Color',[0.5020 0 0.5020]);
text(5.356-1.76, 1.00+0.009, '$(5.356,\;1.00)$', ...
    'Interpreter','latex','FontSize',20, ...
    'Color','red');
text(5.459+0.2, 0.9996-0.009, '$(5.459,\;0.9996)$', ...
    'Interpreter','latex','FontSize',20, ...
    'Color','black');
text(8.929-0.9, 1.0002+0.009, '$(8.929,\;1.0002)$', ...
    'Interpreter','latex','FontSize',20, ...
    'Color','black');
annotation('arrow', [0.6 0.65], [0.754 0.6], ...
    'Color','red','LineWidth',2,'HeadStyle','vback2');
title('\textbf{Zoom View}','Interpreter','latex','FontSize',20);

%-----------------------------------------------------------------------
figure(2); clf
grid on; hold on        
axMain = gca;
hold(axMain,'on'); grid(axMain,'on')
plot(axMain, t_tr, dx2dt_2, ...
    'Color',[1 0.5 0],'LineWidth',4);
title(axMain,'\textbf{Lifting Unit Position}', ...
    'Interpreter','latex','FontSize',30);
ylabel(axMain,'$x_2(t)\; \mathrm{(m)}$', ...
    'Interpreter','latex','FontSize',20);
xlabel(axMain,'Time (s)', ...
    'Interpreter','latex','FontSize',20);
axis(axMain,[0 10 0 0.45]);
plot(axMain, 0.382,0.1251,'o', ...
    'MarkerSize',7,'LineWidth',1.5, ...
    'MarkerEdgeColor',[0.7804, 0.0824, 0.5216]);
text(axMain, 0.382+0.3, 0.1251, ...
    '$(0.382,\;0.1251)$', ...
    'Interpreter','latex','FontSize',20, ...
    'Color',[0.7804 0.0824 0.5216]);

% Zoom inset
zoom_x_start = 2.5;
zoom_x_end   = 10;
zoom_y_start = 0.39;
zoom_y_end   = 0.41;
rectangle(axMain,'Position', ...
    [zoom_x_start, zoom_y_start, ...
     zoom_x_end-zoom_x_start, ...
     zoom_y_end-zoom_y_start], ...
    'EdgeColor','red','LineWidth',2);
axZoom = axes('Position',[0.4 0.2 0.5 0.4]); 
box(axZoom,'on'); grid(axZoom,'on'); hold(axZoom,'on')
plot(axZoom, t_tr, dx2dt_2, ...
    'Color',[1 0.5 0],'LineWidth',3);
xlim(axZoom,[zoom_x_start zoom_x_end]);
ylim(axZoom,[zoom_y_start zoom_y_end]);

set(axZoom,'FontSize',8);
plot(axZoom, 3.082,0.392,'o','MarkerSize',7,'LineWidth',1.5,...
    'MarkerEdgeColor',[0.5020 0 0.5020]);
plot(axZoom, 5.478,0.4,'o','MarkerSize',7,'LineWidth',1.5,...
    'MarkerEdgeColor','red');
plot(axZoom, 5.553,0.3999,'o','MarkerSize',7,'LineWidth',1.5,...
    'MarkerEdgeColor','black');
plot(axZoom, 8.993,0.4001,'o','MarkerSize',7,'LineWidth',1.5,...
    'MarkerEdgeColor','black');

text(axZoom, 3.082+0.2,0.392,'$(3.082,\;0.392)$', ...
    'Interpreter','latex','FontSize',20, ...
    'Color',[0.5020 0 0.5020]);
text(axZoom, 5.478-1.43,0.4+0.0015,'$(5.478,\;0.4)$', ...
    'Interpreter','latex','FontSize',20,'Color','red');
text(axZoom, 5.553+0.2,0.3999-0.0015,'$(5.553,\;0.3999)$', ...
    'Interpreter','latex','FontSize',20,'Color','black');
text(axZoom, 8.993-0.92,0.4001+0.0015,'$(8.993,\;0.4001)$', ...
    'Interpreter','latex','FontSize',20,'Color','black');
annotation('arrow', [0.6 0.65], [0.815 0.6], ...
    'Color','red','LineWidth',2,'HeadStyle','vback2');
title(axZoom,'\textbf{Zoom View}', 'Interpreter','latex','FontSize',20);

%--------------------------------------------------------------------------
figure(3)
grid on;
hold on;
% Độ lắc tại đỉnh thanh
plot(t_tr,w(n,:) - w(1,:),'Color',[0.07,0.62,1.00],'LineWidth',2); 
title('\textbf{Beam''s Top Vibration}', 'Interpreter','latex', 'FontSize', 30);
ylabel('$x_4(t)\ \mathrm{(m)}$', 'Interpreter','latex', 'FontSize', 20);
xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 20);
axis([0 10 -8.5*10^-3 8.5*10^-3]);

% Zoom inset
zoom_x_start = 3;
zoom_x_end   = 4;
zoom_y_start = -7e-3;
zoom_y_end   =  7e-3;
rectangle('Position', ...
    [zoom_x_start, zoom_y_start, ...
     zoom_x_end-zoom_x_start, zoom_y_end-zoom_y_start], ...
    'EdgeColor','red','LineWidth',2);
axes('Position',[0.60 0.47 0.35 0.4]);
box on; hold on; grid on;
plot(t_tr, w(n,:) - w(1,:), 'Color',[0.07,0.62,1.00], 'LineWidth',1.2);
xlim([zoom_x_start zoom_x_end]);
ylim([zoom_y_start zoom_y_end]);
set(gca,'FontSize',8);
annotation('arrow', [0.44 0.6], [0.5 0.68], ...
    'Color','red','LineWidth',2,'HeadStyle','vback2');

y_max =  5.5e-3;
y_min = -5.2e-3;
plot([3 4], [y_max y_max], 'r--', 'LineWidth',2);
plot([3 4], [y_min y_min], 'r--', 'LineWidth',2);
x_arrow = 3.5;
line([x_arrow x_arrow], [y_min y_max], 'Color','r','LineWidth',1.5);
plot(x_arrow, 5.3e-3, '^r', 'MarkerFaceColor','r');
plot(x_arrow, -5e-3, 'vr', 'MarkerFaceColor','r');
text(x_arrow-0.15, 0, '$10.6\times10^{-3}\ \mathrm{m}$', ...
    'Interpreter','latex','FontSize',20, ...
    'BackgroundColor','white','EdgeColor','red');
title('\textbf{Zoom View}','Interpreter','latex','FontSize',20);

%--------------------------------------------------------------------------
figure(4); clf
grid on; hold on
axMain = gca;
hold(axMain,'on'); grid(axMain,'on')
plot(axMain, t_tr, F1(1:r), ...
    'Color',[0.07,0.62,1.00],'LineWidth',2.5);
title(axMain,'\textbf{Driving Unit''s Actuated Force}', ...
    'Interpreter','latex','FontSize',30);
ylabel(axMain,'$F_1(t)\ \mathrm{(N)}$', ...
    'Interpreter','latex','FontSize',20);
xlabel(axMain,'Time (s)', ...
    'Interpreter','latex','FontSize',20);
axis(axMain,[0 10 -10 65]);
plot(axMain,0.0040004,59.9599,'o','MarkerSize',10,...
    'MarkerEdgeColor','red','LineWidth',3);
plot(axMain,1.070,-8.0297,'o','MarkerSize',10,...
    'MarkerEdgeColor',[0.5020 0 0.5020],'LineWidth',3);
text(axMain, 0.0040004+0.2, 59.9599, ...
    '$(0.004,\;59.9599)$', ...
    'Interpreter','latex', ...
    'FontSize',20, ...
    'Color','red');
text(axMain, 1.070-0.55, -8.0297+8.2, ...
    '$(1.070,\;-8.0297)$', ...
    'Interpreter','latex', ...
    'FontSize',20, ...
    'Color',[0.5020 0 0.5020]);

% Zoom inset
zoom_x_start = 5;
zoom_x_end   = 6;
zoom_y_start = -1;
zoom_y_end   =  1;
rectangle('Position', ...
    [zoom_x_start, zoom_y_start, ...
     zoom_x_end-zoom_x_start, ...
     zoom_y_end-zoom_y_start], ...
    'EdgeColor',[1 0.5 0], 'LineWidth',2);
axZoom = axes('Position',[0.5 0.45 0.35 0.4]);
box on; hold on; grid on;
plot(t_tr, F1(1:r), 'Color',[0.07 0.62 1.00], 'LineWidth',1.6);
xlim([zoom_x_start zoom_x_end]);
ylim([zoom_y_start zoom_y_end]);
set(gca,'FontSize',8);
title('\textbf{Zoom View}', ...
    'Interpreter','latex','FontSize',14);
y_max =  0.55;
y_min = -0.57;
plot([zoom_x_start zoom_x_end], [y_max y_max], ...
    '--','Color',[1 0.5 0],'LineWidth',2);
plot([zoom_x_start zoom_x_end], [y_min y_min], ...
    '--','Color',[1 0.5 0],'LineWidth',2);
x_arrow = 5.5;
line([x_arrow x_arrow], [y_min y_max], ...
    'Color',[1 0.5 0],'LineWidth',1.8);
plot(x_arrow, y_max-0.05, '^', ...
    'Color',[1 0.5 0],'MarkerFaceColor',[1 0.5 0]);
plot(x_arrow, y_min+0.05, 'v', ...
    'Color',[1 0.5 0],'MarkerFaceColor',[1 0.5 0]);
text(x_arrow-0.05, 0, ...
    '$1.1\ \mathrm{N}$', ...
    'Interpreter','latex','FontSize',20, ...
    'BackgroundColor','white', ...
    'EdgeColor',[1 0.5 0], ...
    'LineWidth',1.5);
annotation('arrow', ...
    [0.555 0.675], [0.23 0.437], ...
    'Color',[1 0.5 0],'LineWidth',2,'HeadStyle','vback2');

%--------------------------------------------------------------------------
% figure(5)
% plot(t_tr,x_IS(:));

%--------------------------------------------------------------------------
figure(5)
Fs = 1/delta_t; 
time = 0:delta_t:(r - 2)*delta_t - delta_t; 
l1 = r; 
fft_w = fft(w(n,:) - w(1,:),l1)*(2/l1); 
abs_w = abs(fft_w); 
freq = 0:(1/time(end)):Fs/2 - (1/time(end)); 
plot(freq,abs_w(1:length(freq)),'Color',[0.07,0.62,1.00],'LineWidth',2); 
hold on; grid on; 
title('\textbf{Model''s Natural Frequency}', ... 
    'Interpreter','latex','FontSize',30); 
xlabel('$f_n\ \mathrm{(Hz)}$', ... 
    'Interpreter','latex','FontSize',20); 
ylabel('Amplitude', ... 
    'Interpreter','latex','FontSize',20); 
axis([0 16 0 5.5*10^-3])

axFFT = gca;  
hold(axFFT,'on');
f1 = 0.3334;    A1 = 0.000153187;
f2 = 5.60112;   A2 = 0.00512493;
f3 = 14.1362;   A3 = 0.000258639;
plot(axFFT,f1,A1,'ko','MarkerSize',6,'LineWidth',1.5);
plot(axFFT,f2,A2,'ro','MarkerSize',6,'LineWidth',1.8);
plot(axFFT,f3,A3,'ko','MarkerSize',6,'LineWidth',1.5);
text(axFFT, f1+0.4, A1+0.00025, ...
    '$f = 0.3334\ \mathrm{Hz}$', ...
    'Interpreter','latex', ...
    'FontSize',18, ...
    'Color','black', ...
    'BackgroundColor','white', ...
    'EdgeColor','black', ...
    'LineWidth',1.2);
text(axFFT, f2+0.5, A2, ...
    '$f = 5.60112\ \mathrm{Hz}$', ...
    'Interpreter','latex', ...
    'FontSize',20, ...
    'Color','red', ...
    'BackgroundColor','white', ...
    'EdgeColor','red', ...
    'LineWidth',1.6);
text(axFFT, f3-2.5, A1+0.00025, ...
    '$f = 14.1362\ \mathrm{Hz}$', ...
    'Interpreter','latex', ...
    'FontSize',18, ...
    'Color','black', ...
    'BackgroundColor','white', ...
    'EdgeColor','black', ...
    'LineWidth',1.2);
