clc;
close all;
clear;
% Thông số của dầm
pA = 0.297; EI = 0.754; L = 0.63; 

% Các vật nặng
mw = 13.1; mk = 0.04;

% Thông số không gian và thời gian
n = 7; r = 500000;
tmax = 15;
dt = tmax/(r - 1);
ds = L/(n - 1); 

% Lực tác động vào xe con
F1 = zeros(1,r);
F1(1:r/2) = 10;

% Lực tác động vào xe nâng
F2 = zeros(1,r);
F2(1:r/2) = 10;

w = zeros(n,r);
%--------------------------------------------------------------------------
for j = 2:(r - 1)
    wsss0 = (w(3,j) - 2*w(2,j) + w(1,j))/(2*ds^3);
    w(1,j + 1) = (F1(j + 1) - EI*wsss0)*(dt^2/mw) + 2*w(1,j) - w(1,j - 1);
    for i = 3:(n - 2)
        S2 = (-EI*dt^2)/(ds^4*pA);
        wssss = w(i + 2,j) - 4*w(i + 1,j) + 6*w(i,j) - 4*w(i - 1,j) + w(i - 2,j);
        w(i,j + 1) = S2*wssss + 2*w(i,j) - w(i,j - 1);
    end
    S3 = (EI*(dt^2))/(mk*2*ds^3);
    w(2,j + 1) = w(1,j + 1);
    wsssl = (-2*w(n,j) + 3*w(n - 1,j) - w(n - 2,j));
    w(n,j + 1) = 2*w(n,j) - w(n,j - 1) + S3*wsssl;
    w(n - 1,j + 1) = (w(n,j + 1) + w(n - 2,j + 1))/2;
end
%--------------------------------------------------------------------------
w0 = linspace(0,1,n);
t_tr = linspace(0,tmax,r);

figure(1)
grid on;
hold on;
plot(t_tr,w(1,:));
title({'Vị trí của xe con'});
ylabel('x1','FontSize',12);
xlabel('t','FontSize',12);

figure(2)
grid on;
hold on;
plot(t_tr,w(n,:) - w(1,:)); % Vị trí tương đối của thanh so với xe con
title({'Vị trí của đỉnh thanh'});
ylabel('x4','FontSize',12);
xlabel('t','FontSize',12);