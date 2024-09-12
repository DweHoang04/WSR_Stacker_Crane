clc; clear;

% Thiết lập các thông số của mô hình
rho_A = 2.1; EI = 14.97;
L = 1.54; l = 2;
mw = 13.1; mh = 0.87; mk = 0.19; m = 14.15;
tmax = 15;
n = 30;
r = 30000;
g = 9.8;
sp2 = 1;
e = zeros(1,r);
S = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                              

dt = tmax / (r - 1); % Độ lớn phần tử thời gian
ds = L / (n - 1); % Độ dài phần tử thanh

S1 = (EI * (dt^2)) / (mw * (ds^3));
S2 = (-EI * (dt^2)) / ((ds^4) * rho_A);
S3 = (EI * (dt^2)) / (mk * 2 * ds^3);

% Input của mô hình
F1 = zeros(1,r); % Lực tác động vào xe con
F2 = zeros(1,r); % Lực tác động vào xe nâng
x2 = zeros(1,r); % Tọa độ của xe nâng theo phương Oy
x3 = zeros(1,r); % Tọa độ của xe nâng theo phương Ox
F1(1:r/2) = 10; % Chưa rõ dùng để làm gì?
F1(r/2:r) = 10;
F2(1:r/2) = 8; F2(r/2:r) = 0;
p = zeros(1,r); 
w = zeros(n,r);

% Mô phỏng điểm x1
for j = 2:(r - 1)
wsss0 = (w(3,j) - 2 * w(2,j) + w(1,j)) / (2 * ds^3); 
w(1,j + 1) = (F1(j + 1) - EI * wsss0) * (dt^2 / mw) + 2 * w(1,j) - w(1,j - 1);
    for i = 3:n - 2
        S2 = (-EI * dt^2) / (ds^4 * rho_A);
        wssss = w(i + 2,j) - 4 * w(i + 1,j) + 6 * w(i,j) - 4 * w(i - 1,j) + w(i - 2,j);
        w(i,j + 1) = S2 * wssss + 2 * w(i,j) - w(i,j - 1);
    end
    w(2,j + 1) = w(1,j + 1);
    wsssl = (-2 * w(n,j) + 3 * w(n - 1,j) - w(n - 2,j));
    w(n,j + 1) = 2 * w(n,j) - w(n,j - 1) + S3 * wsssl; %Dao động điểm cuối thanh
    w(n - 1,j + 1) = (w(n,j + 1) + w(n - 2,j + 1)) / 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                              

    % PID
    e(j) = sp2 - x2(j); % Sai lệch 
    % Thiết lập thông số PID
    kp = 10;
    ki = 8;
    kd = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Điều khiển bằng PID
    S = S + e(j) * dt;
    F2(j + 1) = kp * e(j) + ki * S + kd * (e(j) - e(j - 1))/dt;
    if F2(j + 1) > 15
        F2(j + 1) = 15;
    end
    if F2(j + 1) < -15
        F2(j + 1) = -15;
    end
    % Tìm vị trí xe nâng
    ddtx3 = (x3(j + 1) - 2 * x3(j) + x3(j - 1)) / (dt^2);
    dyx2 = (w(l,j) - w(l - 1,j)) / ds;
    dyx2 = 0;
    % if ddtx3*dyx2*mh > 0.001
    %    dyx2=0;
    % end
    %    dyx2=0;
    x2(j + 1) = (F2(j + 1) - mh * g - mh * ddtx3 * dyx2) * dt^2 / mh + 2 * x2(j) - x2(j - 1);
    l = ceil(x2(j + 1) / ds);
    if x2(j + 1) < ds
        l = 2;
        x2(j + 1) = ds;
    end
    if x2(j + 1) > (L-ds)
        l = n - 1;
        x2(j + 1) = L - ds;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                              

    % Độ lắc xe nâng
    if l == 2 
        x3(j + 1) = w(1,j + 1);
    end
    if l > 2 && l < n -1
        x3(j + 1) = -(EI * dt^2) * (w(l + 1,j + 1) - 2 * w(l,j + 1) + w(l - 1,j + 1)) / (3 * ds^3 * mh) + 2 * x3(j) - x3(j - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                              
 
        if abs(x3(j + 1) - w(1,j + 1)) > 1.54 && x3(j + 1) < w(1,j + 1)
            x3(j + 1) = w(1,j + 1) - 1.54; 
        end 
        if abs(x3(j + 1) - w(1,j + 1)) > 1.54 && x3(j + 1) > w(1,j + 1)
            x3(j + 1) = w(1,j + 1) + 1.54;
        end
        % Với p(j) = l;
        if abs(w(l,j + 1) - x3(j + 1)) < 0.012 
            w(l,j + 1) = x3(j + 1);
        end
    end
    if l == (n - 1)
        x3(j + 1) = w(n - 1,j + 1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                              

% Đồ thị chống rung sử dụng phương pháp PID
x = 0:dt:tmax;
figure(1)
grid on 
hold on
plot(x,w(1,:),'b');
title('Vị trí xe con');
figure(3)
hold on 
grid on
plot(x,w(n,:)-w(1,:),'b');
title('Độ lắc điểm cuối');
figure(2)
plot(x,x2(:));
title('Vị trí xe nâng');
check = x3(:) - w(1,:)';
figure(4)
plot(x,check,'b');
title('Độ lắc xe nâng');
% x3(1000:1005)
% w(4,1000:1005)
% x3(1000) - w(1,1000)
% figure(1)
% grid on 
% hold on
% plot(x,w(1,:),'k');
% title('Vi tri xe con');
% figure(3)
% hold on 
% grid on
% plot(x,w(n,:) - w(1,:),'k');
% title('Diem cuoi');
% figure(2)
% plot(x,x2(:),'k');
% check=x3(:) - w(1,:)';
% figure(4)
% plot(x,check,'k');
% title('xe nang');