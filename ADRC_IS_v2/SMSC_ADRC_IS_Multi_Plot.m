clc;
close all;
clear;

%--------------------------------------------------------------------------
% 1. CẤU HÌNH CÁC TRƯỜNG HỢP MUỐN SO SÁNH
%--------------------------------------------------------------------------
% Cấu trúc: {Tên hiển thị, màu sắc, flag, flag_is}
% flag = 1: ADRC; flag = 2: ADRC + IS
% flag_is = 1: ZV; flag_is = 2: ZVD; flag = 3: ETM4
Scenarios = {
    'ADRC',             [1 0.5 0], 1, 0; % Light Blue [0.07,0.62,1.00] Orange [1 0.5 0]
    'ADRC + ZV',        [0 0.45 0.74], 2, 1; % Red/Orange [0.8500 0.3250 0.0980] Dark Blue [0 0.45 0.74]
    % 'ADRC + ZVD',       [0.9290 0.6940 0.1250], 2, 2; % Yellow
    % 'ADRC + ETM4',      [0.4940 0.1840 0.5560], 2, 3; % Purple
};

num_cases = size(Scenarios, 1);
Results = struct(); % Biến lưu kết quả của tất cả các lần chạy

%--------------------------------------------------------------------------
% 2. VÒNG LẶP MÔ PHỎNG (CHẠY TỪNG BỘ ĐIỀU KHIỂN)
%--------------------------------------------------------------------------

for k = 1:num_cases
    % Lấy tham số cho lần chạy này
    CaseName = Scenarios{k, 1};
    flag = Scenarios{k, 3};
    flag_is = Scenarios{k, 4};
    
    fprintf('Dang chay mo phong truong hop: %s ...\n', CaseName);

    % --- KHỞI TẠO THÔNG SỐ ---
    L = 1; EI = 14.97; rho_A = 2.1;
    mw = 15; mk = 0.2; mh = 0.9; g = 9.81;
    n = 9; r = 15000;
    tmax = 15;
    delta_Y = L/(n - 1);
    delta_t = tmax/(r - 1);

    w = zeros(n,r);
    x3 = zeros(1,r); 
    x2 = 2;
    dx2dt_2 = zeros(1,r);
    dx2dt_2(1:3) = delta_Y;

    F1 = zeros(1,r); F1(1:r) = 10;
    F2 = zeros(1,r);

    % Thông số ADRC
    T_set = 3; T_sample = delta_t; s_CL = -6/T_set;
    Kp = s_CL^2; Kd = -2*s_CL; K_ESO = 10; b01 = 1/mw; b02 = 1/mh;

    s_ESO = K_ESO*s_CL; z_ESO = exp(s_ESO*T_sample);
    l1 = 1 - z_ESO^3; 
    l2 = (3*(1 - z_ESO)^2*(1 + z_ESO))/(2*T_sample);
    l3 = (1 - z_ESO)^3/(T_sample^2);

    Ad = [1 T_sample (T_sample^2)/2; 0 1 T_sample; 0 0 1];
    Bd1 = [b01*(T_sample^2)/2; b01*T_sample; 0];
    Bd2 = [b02*(T_sample^2)/2; b02*T_sample; 0];
    Cd = [1 0 0]; Lc = [l1; l2; l3];

    A_ESO = Ad - Lc*Cd*Ad;
    B_ESO1 = Bd1 - Lc*Cd*Bd1;
    B_ESO2 = Bd2 - Lc*Cd*Bd2;

    xd = zeros(3,r); xl = zeros(3,r);
    x1_set = 1; x2_set = 0.4;

    % Thông số bộ tạo dạng tín hiệu
    f = 5.60112; 
    zeta = 0; K = exp((zeta*pi)/sqrt(1 - zeta^2)); 
    mopt = 1; 
    t2 = 1/(2*f); t3 = 1/f;
    x_IS = zeros(1,r);

    % --- BẮT ĐẦU VÒNG LẶP TÍNH TOÁN ---
    for j = 3:(r - 1)
        wyyy0 = (w(3,j) - 2*w(2,j) + w(1,j))/(2*delta_Y^3);
        S2 = (F1(j + 1) - EI*wyyy0)/mw;
        w(1,j + 1) = 2*w(1,j) - w(1,j - 1) + delta_t^2*S2; 
        for i = 3:(n - 2)
            wyyyy = w(i + 2,j) - 4*w(i + 1,j) + 6*w(i,j) - 4*w(i - 1,j) + w(i - 2,j);
            S1 = (-EI/(rho_A*delta_Y^4))*wyyyy;
            S4 = (-EI/(mh*delta_Y^3))*wyyyy;
            
            dx3dt_2 = (x3(j) - 2*x3(j - 1) + x3(j - 2))/(2*delta_t^2);
            wyx2 = (w(x2 + 1,j - 1) - w(x2 - 1,j - 1))/(2*delta_Y^2);
            dx2dt_2(j + 1) = 2*dx2dt_2(j) - dx2dt_2(j - 1) + ((F2(j + 1) - mh*g ...
                           - mh*dx3dt_2*wyx2)*delta_t^2)/mh; 

            if dx2dt_2(j + 1) < delta_Y
                x2 = 2; dx2dt_2(j + 1) = delta_Y;
            elseif dx2dt_2(j + 1) > L - delta_Y
                x2 = n - 1; dx2dt_2(j + 1) = L - delta_Y;
            else
                x2 = ceil(dx2dt_2(j + 1)/delta_Y);
            end

            if dx2dt_2(j + 1) ~= x2
                w(i,j + 1) = 2*w(i,j) - w(i,j - 1) + delta_t^2*S1; 
                x3(j + 1) = w(i,j + 1) - w(1,j + 1);
            else
                w(i,j + 1) = 2*w(i,j) - w(i,j - 1) + delta_t^2*S4; 
                x3(j + 1) = w(i,j + 1) - w(1,j + 1);
            end
        end
        wyyyl = (-2*w(n,j) + 3*w(n - 1,j) - w(n - 2,j))/(2*delta_Y^3);
        S3 = (EI/mk)*wyyyl;
        w(2,j + 1) = w(1,j + 1);
        w(n,j + 1) = 2*w(n,j) - w(n,j - 1) + delta_t^2*S3; 
        w(n - 1,j + 1) = (w(n,j + 1) + w(n - 2,j + 1))/2;

        % --- LOGIC ĐIỀU KHIỂN ---
        if flag == 1 % Chỉ ADRC
             % Xe con    
            xd(:,j + 2) = A_ESO*xd(:,j + 1) + B_ESO1*F1(j + 1) + Lc*w(1,j + 1);
            F1(j + 2) = (Kp*(x1_set - xd(1,j + 2)) - Kd*xd(2,j + 2) - xd(3,j + 2))/b01;
            % Xe nâng
            xl(:,j + 2) = A_ESO*xl(:,j + 1) + B_ESO2*(F2(j + 1) - mh*g) + Lc*dx2dt_2(j + 1);
            F2(j + 2) = (Kp*(x2_set - xl(1,j + 2)) - Kd*xl(2,j + 2) - xl(3,j + 2))/b02 + mh*g;
        elseif flag == 2 % ADRC + IS
             if flag_is == 1 % ZV
                A1 = K/(1 + K); A2 = A1;
                x1_set_is = x1_set*(A1*unit_step(j*delta_t) + A2*unit_step(j*delta_t - t2));
            elseif flag_is == 2 % ZVD
                A1 = (K/(1 + K))^2; A2 = (2*K)/(1 + K)^2; A3 = 1/(1 + K)^2;
                x1_set_is = x1_set*(A1*unit_step(j*delta_t) + A2*unit_step(j*delta_t - t2) ...
                                  + A3*unit_step(j*delta_t - t3));
            elseif flag_is == 3 % ETM4
                t2_etm = 1/(3*f); t3_etm = 2/(3*f); t4_etm = 1/f;
                I = (K^2*(1 + mopt))/(K^2 + (1 + mopt)*(K^(4/3) + K^(2/3)) + mopt);
                A1 = I/(1 + mopt); A2 = I/K^(2/3); A3 = I/K^(4/3); A4 = (mopt*I)/(K^2*(1 + mopt));
                x1_set_is = x1_set*(A1*unit_step(j*delta_t) + A2*unit_step(j*delta_t - t2_etm) ...
                                  + A3*unit_step(j*delta_t - t3_etm)) + A4*unit_step(j*delta_t - t4_etm);
            end
            x_IS(j) = x1_set_is;
            
            % Xe con
            xd(:,j + 2) = A_ESO*xd(:,j + 1) + B_ESO1*F1(j + 1) + Lc*w(1,j + 1);
            F1(j + 2) = (Kp*(x1_set_is - xd(1,j + 2)) - Kd*xd(2,j + 2) - xd(3,j + 2))/b01;
            % Xe nâng
            xl(:,j + 2) = A_ESO*xl(:,j + 1) + B_ESO2*(F2(j + 1) - mh*g) + Lc*dx2dt_2(j + 1);
            F2(j + 2) = (Kp*(x2_set - xl(1,j + 2)) - Kd*xl(2,j + 2) - xl(3,j + 2))/b02 + mh*g;
        end
    end
    
    % --- LƯU KẾT QUẢ VÀO STRUCT ---
    Results(k).w = w;
    Results(k).dx2dt_2 = dx2dt_2;
    Results(k).top_vibration = w(n,:) - w(1,:);
    Results(k).F1 = F1(1:r);
end

%--------------------------------------------------------------------------
% 3. VẼ ĐỒ THỊ SO SÁNH
%--------------------------------------------------------------------------
t_tr = 0:delta_t:tmax;

% --- Figure 1: Driving Unit Position ---
figure(1); clf; hold on; grid on;
% Vị trí xe con của hai bộ điều khiển
for k = 1:num_cases
    plot(t_tr, Results(k).w(1,:), 'Color', Scenarios{k,2}, 'LineWidth', 4);
end
title('\textbf{Driving Unit Position}', 'Interpreter','latex', 'FontSize', 30);
ylabel('$x_1(t)\ \mathrm{(m)}$', 'Interpreter','latex', 'FontSize', 20);
xlabel('Time (s)','Interpreter','latex','FontSize', 20);
legend(Scenarios(:,1), 'FontSize', 16, 'Interpreter', 'latex', 'Location', 'best');
axis([0 10 0 1.2]);

% Zoom inset
zoom_x_start = 2; zoom_x_end = 4;
zoom_y_start = 0.95; zoom_y_end = 1;
rectangle('Position', ...
    [zoom_x_start, zoom_y_start, ...
     zoom_x_end-zoom_x_start, zoom_y_end-zoom_y_start], ...
    'EdgeColor','red','LineWidth',2);
axes('Position',[0.3 0.2 0.4 0.3]);
box on; hold on; grid on;
for k = 1:num_cases
    plot(t_tr, Results(k).w(1,:), 'Color', Scenarios{k,2}, 'LineWidth', 3);
end
xlim([zoom_x_start zoom_x_end]);
ylim([zoom_y_start zoom_y_end]);
set(gca,'FontSize',8);
p1 = plot(2.915,0.98,'o','MarkerSize',7,'MarkerEdgeColor',[0.07,0.62,1.00],'LineWidth',3,'MarkerIndices',1);
p2 = plot(2.959,0.98,'o','MarkerSize',7,'MarkerEdgeColor',[0.8500 0.3250 0.0980],'LineWidth',3,'MarkerIndices',1);
text(2.915-0.45, 0.98+0.005, '$(2.915,\;0.98)$', ...
    'Interpreter','latex','FontSize',20, ...
    'Color',[0.07,0.62,1.00]);
text(2.959+0.07, 0.98-0.005, '$(2.959,\;0.98)$', ...
    'Interpreter','latex','FontSize',20, ...
    'Color',[0.8500 0.3250 0.0980]);
annotation('arrow', [0.363 0.5], [0.754 0.5], ...
    'Color','red','LineWidth',2,'HeadStyle','vback2');
title('\textbf{Zoom View}','Interpreter','latex','FontSize',20);

% --- Figure 2: Lifting Unit Position ---
figure(2); clf; hold on; grid on;
% Vị trí xe nâng của hai bộ điều khiển
for k = 1:num_cases
    plot(t_tr, Results(k).dx2dt_2, 'Color', Scenarios{k,2}, 'LineWidth', 4);
end
title('\textbf{Lifting Unit Position}', 'Interpreter','latex', 'FontSize', 30);
ylabel('$x_2(t)\ \mathrm{(m)}$', 'Interpreter','latex', 'FontSize', 20);
xlabel('Time (s)','Interpreter','latex','FontSize', 20);
legend(Scenarios(:,1), 'FontSize', 16, 'Interpreter', 'latex', 'Location', 'best');
axis([0 10 0 0.45]);

% Zoom inset
zoom_x_start = 2.5; zoom_x_end = 4.5;
zoom_y_start = 0.385; zoom_y_end = 0.4;
rectangle('Position', ...
    [zoom_x_start, zoom_y_start, ...
     zoom_x_end-zoom_x_start, zoom_y_end-zoom_y_start], ...
    'EdgeColor','red','LineWidth',2);
axes('Position',[0.3 0.2 0.4 0.3]);
box on; hold on; grid on;
for k = 1:num_cases
    plot(t_tr, Results(k).dx2dt_2, 'Color', Scenarios{k,2}, 'LineWidth', 3);
end
xlim([zoom_x_start zoom_x_end]);
ylim([zoom_y_start zoom_y_end]);
set(gca,'FontSize',8);
p11 = plot(3.082,0.392,'o','MarkerSize',7,'MarkerEdgeColor',[1 0.5 0],'LineWidth',3,'MarkerIndices',1);
p21 = plot(3.166,0.392,'o','MarkerSize',7,'MarkerEdgeColor',[0 0.45 0.74],'LineWidth',3,'MarkerIndices',1);
text(3.082-0.45, 0.392+0.001, '$(3.082,\;0.392)$', ...
    'Interpreter','latex','FontSize',20, ...
    'Color',[1 0.5 0]);
text(3.166+0.07, 0.392-0.001, '$(3.166,\;0.392)$', ...
    'Interpreter','latex','FontSize',20, ...
    'Color',[0 0.45 0.74]);
annotation('arrow', [0.4 0.5], [0.808 0.5], ...
    'Color','red','LineWidth',2,'HeadStyle','vback2');
title('\textbf{Zoom View}','Interpreter','latex','FontSize',20);

% --- Figure 3: Beam Top Vibration ---
figure(3); clf; hold on; grid on;
% Độ lắc tại đỉnh thanh
for k = 1:num_cases
    plot(t_tr, Results(k).top_vibration, 'Color', Scenarios{k,2}, 'LineWidth', 2);
end
title('\textbf{Beam''s Top Vibration}', 'Interpreter','latex', 'FontSize', 30);
ylabel('$x_4(t)\ \mathrm{(m)}$', 'Interpreter','latex', 'FontSize', 20);
xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 20);
legend(Scenarios(:,1), 'FontSize', 16, 'Interpreter', 'latex', 'Location', 'best');
y_max =  5.4e-3;
plot([0 10], [y_max y_max], 'r--', 'LineWidth', 2, 'HandleVisibility', 'off');
x_arrow = 3.25;
line([x_arrow x_arrow], [2e-3 y_max], 'Color','r','LineWidth',1.5, 'HandleVisibility','off');
plot(x_arrow, 5.3e-3, '^r', 'MarkerFaceColor','r', 'HandleVisibility','off');
text(x_arrow-0.5, 2e-3, '$5.2\times10^{-3}\ \mathrm{m}$', ...
    'Interpreter','latex','FontSize',20, ...
    'BackgroundColor','white','EdgeColor','red');
p12 = plot(0.0860057,-0.00402679,'o','MarkerSize',7,'MarkerEdgeColor','black', ...
        'LineWidth',3,'MarkerIndices',1,'HandleVisibility','off');
p22 = plot(0.0830055,-0.00810709,'o','MarkerSize',7,'MarkerEdgeColor','black', ...
        'LineWidth',3,'MarkerIndices',1,'HandleVisibility','off');
axis([0 10 -8.5*10^-3 8.5*10^-3]);

% Zoom inset
zoom_x_start = 3;
zoom_x_end   = 3.5;
zoom_y_start = -4.5e-4;
zoom_y_end   =  5.5e-4;
rectangle('Position', ...
    [zoom_x_start, zoom_y_start, ...
     zoom_x_end-zoom_x_start, zoom_y_end-zoom_y_start], ...
    'EdgeColor','black','LineWidth',2);
axes('Position',[0.57 0.2 0.35 0.4]);
box on; hold on; grid on;
for k = 1:num_cases
    plot(t_tr, Results(k).top_vibration, 'Color', Scenarios{k,2}, 'LineWidth', 1.2);
end
xlim([zoom_x_start zoom_x_end]);
ylim([zoom_y_start zoom_y_end]);
set(gca,'FontSize',8);
annotation('arrow', [0.402 0.57], [0.52 0.38], ...
    'Color','black','LineWidth',2,'HeadStyle','vback2');

y_max =  4.6e-4;
plot([3 4], [y_max y_max], 'black--', 'LineWidth',2);
x_arrow = 3.25;
line([x_arrow x_arrow], [0 y_max], 'Color','black','LineWidth',1.5);
plot(x_arrow, 4.45e-4, '^black', 'MarkerFaceColor','black');
text(x_arrow-0.055, 0, '$4.6\times10^{-4}\ \mathrm{m}$', ...
    'Interpreter','latex','FontSize',20, ...
    'BackgroundColor','white','EdgeColor','black');
title('\textbf{Zoom View}','Interpreter','latex','FontSize',20);

% --- Figure 4: Driving Force ---
figure(4); clf; ax_main = axes;
hold(ax_main, 'on'); grid(ax_main, 'on');
for k = 1:num_cases
    plot(ax_main, t_tr, Results(k).F1, 'Color', Scenarios{k,2}, 'LineWidth', 4);
end
title(ax_main, '\textbf{Driving Unit''s Actuated Force}', 'Interpreter','latex','FontSize',20);
ylabel(ax_main, '$F_1(t)\ \mathrm{(N)}$', 'Interpreter','latex','FontSize',20);
xlabel(ax_main, 'Time (s)', 'Interpreter','latex','FontSize',20);
legend(ax_main, Scenarios(:,1), 'FontSize', 16, 'Interpreter', 'latex');
axis(ax_main, [0 10 -10 65]);

% --- Zoom inset 1 ---
zoom_x_start = 0;
zoom_x_end   = 1;
zoom_y_start = -10;
zoom_y_end   =  65;

rectangle(ax_main, 'Position', ... 
    [zoom_x_start, zoom_y_start, ...
     zoom_x_end-zoom_x_start, zoom_y_end-zoom_y_start], ...
    'EdgeColor','black','LineWidth',2);
axZoom1 = axes('Position',[0.3 0.6 0.25 0.27]);
box on; hold on; grid on;
for k = 1:num_cases
    plot(t_tr, Results(k).F1, 'Color', Scenarios{k,2}, 'LineWidth', 3);
end
xlim([zoom_x_start zoom_x_end]);
ylim([zoom_y_start zoom_y_end]);
set(gca,'FontSize',8);
annotation('arrow', [0.208 0.3], [0.53 0.63], ...
    'Color','black','LineWidth',2,'HeadStyle','vback2');
plot(axZoom1, 0.00400027,59.9599,'o','MarkerSize',7,'LineWidth',3,...
    'MarkerEdgeColor',[0.07,0.62,1.00]);
text(axZoom1, 0.00400027+0.05,59.9599-3,'$(0.004,\;59.9599)$', ...
    'Interpreter','latex','FontSize',20, ...
    'Color',[0.07,0.62,1.00]);
plot(axZoom1, 0.0910061,50.3158,'o','MarkerSize',7,'LineWidth',3,...
    'MarkerEdgeColor',[0.8500 0.3250 0.0980]);
text(axZoom1, 0.0910061+0.05,50.3158-3,'$(0.091,\;50.3158)$', ...
    'Interpreter','latex','FontSize',20, ...
    'Color',[0.8500 0.3250 0.0980]);
title('\textbf{Zoom View 1}','Interpreter','latex','FontSize',20);

% --- Zoom inset 2 ---
zoom_x_start1 = 5;
zoom_x_end1   = 6;
zoom_y_start1 = -1;
zoom_y_end1   =  1;

rectangle(ax_main, 'Position', ...  
    [zoom_x_start1, zoom_y_start1, ...
     zoom_x_end1-zoom_x_start1, zoom_y_end1-zoom_y_start1], ...
    'EdgeColor','black','LineWidth',2);
axZoom2 = axes('Position',[0.65 0.3 0.25 0.27]); % Tạo trục mới
box on; hold on; grid on;
for k = 1:num_cases
    plot(t_tr, Results(k).F1, 'Color', Scenarios{k,2}, 'LineWidth', 3);
end
xlim([zoom_x_start1 zoom_x_end1]);
ylim([zoom_y_start1 zoom_y_end1]);
set(gca,'FontSize',8);
annotation('arrow', [0.555 0.65], [0.23 0.43], ... 
    'Color','black','LineWidth',2,'HeadStyle','vback2');
title('\textbf{Zoom View 2}','Interpreter','latex','FontSize',20);

% --- Helper Functions ---
function y = unit_step(t)
    y = double(t >= 0);
end