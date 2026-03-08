clc; clear; close all;

%% Parameters
% Tham số thang nâng
pA = 1;
EI = 1;
g  = 9.81;

L  = 0.63;
mw = 50;
mh = 7.5;
mk = 0.04;

c = 5; % Hệ số giảm chấn

%% Simulation settings
tmax = 15;
n    = 9;
r    = 15000;

dt = tmax / (r - 1);
ds = L / (n - 1);

% Setpoints
sp1 = 1;
sp2 = 8 * ds + 0.02;

%% State / signal initialization
% Signals
e  = zeros(1, r);
F1 = zeros(1, r);
F2 = zeros(1, r);
h  = zeros(1, r);
w3 = zeros(1, r);

% Initial conditions
h(1:3)   = ds;
F2(1:r)  = mh * g - 0.1;

% l: initial lift-car index on the beam
l = 2;

% Beam states
w  = zeros(n, r);
dw = zeros(1, r);
z1 = zeros(1, r);
dm = zeros(1, r);

%% Controller selection
% flag = 0: open-loop test (forces)
% flag = 1: ROBUST-PD
% flag = 2: barrier
% flag = 3: PD
% flag = 4: disturbance observer
% flag = 5..7: other controllers (kept as-is)
flag = 0;

%% Simulation
for j = 3 : (r - 1)
    % 3rd spatial derivative at the base (as in the original)
    wsss0 = (w(3, j) - 2 * w(2, j) + w(1, j)) / (2 * ds^3);

    % Base mass dynamics
    w(1, j + 1) = (F1(j + 1) - EI * wsss0 + dw(j + 1)) * (dt^2 / mw) ...
        + 2 * w(1, j) - w(1, j - 1);

    for i = 3 : (n - 2)
        % Common terms
        S2  = (-EI * dt^2) / (ds^4 * pA);
        S21 = (-EI * dt^2) / (ds^3 * mh);
        C   = -c * dt^2 * (w(i, j) - w(i, j - 1)) / dt;

        wssss  = w(i + 2, j) - 4 * w(i + 1, j) + 6 * w(i, j) - 4 * w(i - 1, j) + w(i - 2, j);
        ddtx3  = (w3(j) - 2 * w3(j - 1) + w3(j - 2)) / (dt^2);
        dyx2   = ((w(l + 1, j - 1) - w(l - 1, j - 1)) / 2) * ds;

        % Lift-car position dynamics
        h(j + 1) = (F2(j + 1) - mh * g - mh * ddtx3 * dyx2) * dt^2 / mh ...
            + 2 * h(j) - h(j - 1);

        % Saturate and update l
        if h(j + 1) < ds
            l      = 2;
            h(j + 1) = ds;
        elseif h(j + 1) > (L - ds)
            l      = n - 1;
            h(j + 1) = L - ds;
        else
            l = ceil(h(j + 1) / ds);
        end

        % Beam update at each node
        if i == l
            w(i, j + 1) = S21 * wssss + C + 2 * w(i, j) - w(i, j - 1);
        else
            w(i, j + 1) = S2 * wssss + C + 2 * w(i, j) - w(i, j - 1);
        end

        % This is how w3 was computed in the original (overwritten each i)
        w3(j + 1) = w(i, j + 1) - w(1, j + 1);

        % Boundary conditions
        wsssl = (-2 * w(n, j) + 3 * w(n - 1, j) - w(n - 2, j));
        S3    = (EI * (dt^2)) / (mk * 2 * ds^3);

        w(2, j + 1)     = w(1, j + 1);
        w(n, j + 1)     = 2 * w(n, j) - w(n, j - 1) + S3 * wsssl;
        w(n - 1, j + 1) = (w(n, j + 1) + w(n - 2, j + 1)) / 2;

        %% Controller
        switch flag
            case 0
                F1(1:2000)  = 10;
                F1(2001:4000) = -10;
                F2(1:2000)  = 10;
                F2(2001:4000) = -10;
                check2 = zeros(1, r); 

            case 1
                % ROBUST-PD
                k1 = 10;
                kx = 20;
                kn = 1.5;

                dx1     = (w(1, j + 1) - w(1, j)) / dt;
                F1(j + 2) = -k1 * (w(1, j + 1) - sp1) - kx * dx1;

                S2_ = (w(l + 1, j + 1) - 2 * w(l, j + 1) + w(l - 1, j + 1)) / (ds^3);
                dx2 = ((w(l + 1, j + 1) - w(l - 1, j + 1)) / 2) * ds;
                dtx2 = (h(j + 1) - h(j)) / dt;

                F2(j + 2) = mh * g + EI * S2_ * dx2 - (h(j + 1) - sp2) - kn * dtx2;

            case 2
                % barrier
                k1   = 20;
                kmax = 0.004;
                k2   = 0.01;
                k3   = 30;
                k0   = 10;
                kc   = 0.001;
                kn   = 1.5;

                z   = w(l, j + 1);
                dwl = (w(1, j + 1) - w(1, j)) / dt;

                S1 = kc / ((kmax^2 - z^2) * 2.3) + k2;
                F1(j + 2) = -k1 * (w(1, j + 1) - sp1) - k3 * dwl - S1 * kmax * k0 * dwl - S1 * kmax;

                S2_  = (w(l + 1, j + 1) - 2 * w(l, j + 1) + w(l - 1, j + 1)) / (ds^3);
                dx2  = ((w(l + 1, j + 1) - w(l - 1, j + 1)) / 2) * ds;
                dtx2 = (h(j + 1) - h(j)) / dt;

                F2(j + 2) = mh * g + EI * S2_ * dx2 - (h(j + 1) - sp2) - kn * dtx2;

            case 3
                % PD
                ki1 = 1;
                kd1 = 10;

                e(j + 1)  = (sp1 - w(1, j + 1));
                de        = (e(j + 1) - e(j)) / dt;
                F1(j + 2) = ki1 * (sp1 - w(1, j + 1)) + kd1 * de;

            case 4
                % disturbance observer
                a1 = 1;
                k1 = 2.5;
                k2 = 25;
                kn = 1.5;

                dwl = (w(1, j + 1) - w(1, j)) / dt;

                dz1      = -a1 * (F1(j + 1) - EI * wsss0) - a1 * (z1(j + 1) + a1 * mw * dwl);
                z1(j + 1) = dz1 * dt + z1(j);

                d1e       = z1(j + 1) + a1 * mw * dwl;
                F1(j + 2) = -k1 * dwl - d1e - k2 * (w(1, j + 1) - sp1);

                S2_  = (w(l + 1, j + 1) - 2 * w(l, j + 1) + w(l - 1, j + 1)) / (ds^3);
                dx2  = (w(l + 1, j + 1) - w(l - 1, j + 1)) / ds;
                dtx2 = (h(j + 1) - h(j)) / dt;

                F2(j + 2) = mh * g + EI * S2_ * dx2 - (h(j + 1) - sp2) - kn * dtx2;

            case 5
                % (kept from original)
                um = 15;
                vm = 60;
                c1 = 10;
                c3 = 13;
                c2 = 40;
                k1 = 0.12;
                kn = 1.5;

                dw1  = (w(1, j + 1) - w(1, j)) / dt;
                wsss = (w(3, j + 1) - 3 * w(2, j + 1) + 2 * w(1, j + 1)) / ds(j + 1)^3; 

                e(j + 1) = (sp1 - w(1, j + 1));
                de       = (e(j + 1) - e(j)) / dt;

                a1(j + 1) = -c1 * dw1 + 1/2 * EI * wsss + 1/2 * k1 * e(j + 1);
                d11       = dw1;
                d21       = F1(j + 1) - a1(j + 1);

                a2(j + 1) = -2 * (1 / mw) * d11 - (1 / mw) * c2 * d21 + (a1(j + 1) - a1(j)) / dt;
                d31(j + 1) = dF1(j + 1) - a2(j + 1); 

                da2     = (a2(j + 1) - a2(j)) / dt;
                U1(j + 2) = (-1 / mw) * c3 * d31(j + 1) - d21 + da2; 

                dv1       = (1 / (cosh(v1(j + 1) / vm))^2) * U1(j + 2);
                v1(j + 2) = dt * dv1 + v1(j + 1);

                dF1(j + 2) = vm * tanh(v1(j + 2) / vm);
                du1        = (1 / (cosh(u1(j + 1) / vm))^2) * dF1(j + 2);

                u1(j + 2) = du1 * dt + u1(j + 1);
                F1(j + 2) = um * tanh(u1(j + 2) / um);

                S2_  = (w(l + 1, j + 1) - 2 * w(l, j + 1) + w(l - 1, j + 1)) / (ds^3);
                dx2  = ((w(l + 1, j + 1) - w(l - 1, j + 1)) / 2) * ds;
                dtx2 = (h(j + 1) - h(j)) / dt;

                F2(j + 2) = mh * g + EI * S2_ * dx2 - (h(j + 1) - sp2) - kn * dtx2;

            case 6
                % (kept from original)
                kmax = 0.01;
                g1   = 0.1;
                k1   = 1;
                c1   = 1;
                k11  = 1;
                k2   = 1;

                wsss0_ = (w(3, j) - 3 * w(2, j) + 2 * w(1, j)); 
                wsssl_ = (-2 * w(n, j) + 3 * w(n - 1, j) - w(n - 2, j)); 

                dw1 = (w(1, j + 1) - w(1, j)) / dt;
                z(j + 1) = w(n, j + 1) - w(1, j + 1); 

                a = log((2 * kmax^2) / (kmax^2 - z(j + 1)^2)) / log(2.71);
                dz = (z(j + 1) - z(j)) / dt;

                s1 = EI * wsss0_ / ds(j)^3; 

                v1(j + 2) = -k1 * dz - s1 - k11 * dz / a - mw * dz / a * z(j + 1) * dz / (kmax^2 - z(j + 1)^2) ...
                    - (k2 * (w(1, j) - sp1)) / a - mw * EI * wsssl_ / (ds(j)^3 * mk(j) * a);

                p1(j + 2)  = (g1 * dz * v1(j + 2) * a - p1(j + 1)) * dt + p1(j + 1);
                Fc1(j + 2) = p1(j + 2) * v1(j + 2);
                F1(j + 2)  = c1 * Fc1(j + 2);

            case 7
                % ADRC + ZVD xe con / xe nâng (kept from original)
                w1_eso(:, j + 2) = A_ESO1 * w1_eso(:, j + 1) + B_ESO1 * F1(j + 1) + Lc1 * w(1, j + 1);
                F1(j + 2) = (Kp1 * (sp1 * sp1_is(j) - w1_eso(1, j + 2)) - Kd1 * w1_eso(2, j + 2) - w1_eso(3, j + 2)) / b01;

                w2_eso(:, j + 2) = A_ESO2 * w2_eso(:, j + 1) + B_ESO2 * (F2(j + 1) - mh * g) + Lc2 * h(j + 1);
                F2(j + 2) = (Kp2 * (sp2 - w2_eso(1, j + 2)) - Kd2 * w2_eso(2, j + 2) - w2_eso(3, j + 2)) / b02 + mh * g;

            otherwise
                % No controller update
        end
    end
end

%% Plotting
x = 0 : dt : tmax;
y = 0 : ds : L; 

figure(1)
subplot(2, 2, 1)
grid on; hold on;
plot(x, w(1, :), 'b', LineWidth=1);
xlabel('Thời gian(s)');
ylabel('Vị trí xe con(m)');
title('Vị trí xe con');
ylim([0 2]);

subplot(2, 2, 2)
hold on; grid on;

dolac = zeros(1, r);
for j = 1 : r
    dolac(j) = w(n, j) - w(1, j);
end

plot(x, dolac, 'b', LineWidth=1);
title('Độ lắc điểm cuối');
xlabel('Thời gian(s)');
ylabel('m');

subplot(2, 2, 3)
plot(x, h, 'r', LineWidth=1);
grid on;
title('Vị trí xe nâng');
xlabel('Thời gian(s)');
ylabel('m');

subplot(2, 2, 4)
grid on;
plot(x, w3, LineWidth=1);
title('Độ lắc xe nâng');
xlabel('Thời gian(s)');
ylabel('m');

figure(2)
subplot(1, 2, 1)
plot(x, F1(1:r), LineWidth=1.5);
title('F1');
grid on;

check = dolac - w3; 

subplot(1, 2, 2)
plot(x, F2(1:r), LineWidth=1);
title('F2');
grid on;

figure(3)
% Phân tích phổ
Fs   = 1 / dt;
time = 0 : dt : (r - 3) * dt;
l1   = r;

freq = 0 : (1 / time(end)) : Fs/2 - (1 / time(end));

for i = 1 : n
    fft_w = fft(w(i, :) - w(1, :), l1) * (2 / l1);
    abs_w = abs(fft_w);

    subplot(ceil(n / 2), 2, i);
    plot(freq, abs_w(1:length(freq)), 'LineWidth', 0.8);
    hold on; grid on;

    set(gca, 'FontName', 'Times New Roman', 'FontSize', 9);
    xlabel('Frequency(Hz)');
    ylabel('Amplitude(m)');

    title(['FFT của thanh tại điểm ', num2str(i), ' khi xe con ở vị trí ', num2str(l)]);
end
clc; clear; close all;

%% Parameters
% Tham số thang nâng
pA = 1;
EI = 1;
g  = 9.81;

L  = 0.63;
mw = 50;
mh = 7.5;
mk = 0.04;

c = 5; % Hệ số giảm chấn

%% Simulation settings
tmax = 15;
n    = 9;
r    = 15000;

dt = tmax / (r - 1);
ds = L / (n - 1);

% Setpoints
sp1 = 1;
sp2 = 8 * ds + 0.02;

%% State / signal initialization
% Signals
e  = zeros(1, r);
F1 = zeros(1, r);
F2 = zeros(1, r);
h  = zeros(1, r);
w3 = zeros(1, r);

% Initial conditions
h(1:3)   = ds;
F2(1:r)  = mh * g - 0.1;

% l: initial lift-car index on the beam
l = 2;

% Beam states
w  = zeros(n, r);
dw = zeros(1, r);
z1 = zeros(1, r);
dm = zeros(1, r);  % kept for compatibility with the original code

%% Controller selection
% flag = 0: open-loop test (forces)
% flag = 1: ROBUST-PD
% flag = 2: barrier
% flag = 3: PD
% flag = 4: disturbance observer
% flag = 5..7: other controllers (kept as-is)
flag = 4;

%% Simulation
for j = 3 : (r - 2)
    % 3rd spatial derivative at the base (as in the original)
    wsss0 = (w(3, j) - 2 * w(2, j) + w(1, j)) / (2 * ds^3);

    % Base mass dynamics
    w(1, j + 1) = (F1(j + 1) - EI * wsss0 + dw(j + 1)) * (dt^2 / mw) ...
        + 2 * w(1, j) - w(1, j - 1);

    for i = 3 : (n - 2)
        % Common terms
        S2  = (-EI * dt^2) / (ds^4 * pA);
        S21 = (-EI * dt^2) / (ds^3 * mh);
        C   = -c * dt^2 * (w(i, j) - w(i, j - 1)) / dt;

        wssss  = w(i + 2, j) - 4 * w(i + 1, j) + 6 * w(i, j) - 4 * w(i - 1, j) + w(i - 2, j);
        ddtx3  = (w3(j) - 2 * w3(j - 1) + w3(j - 2)) / (dt^2);
        dyx2   = ((w(l + 1, j - 1) - w(l - 1, j - 1)) / 2) * ds;

        % Lift-car position dynamics
        h(j + 1) = (F2(j + 1) - mh * g - mh * ddtx3 * dyx2) * dt^2 / mh ...
            + 2 * h(j) - h(j - 1);

        % Saturate and update l
        if h(j + 1) < ds
            l      = 2;
            h(j + 1) = ds;
        elseif h(j + 1) > (L - ds)
            l      = n - 1;
            h(j + 1) = L - ds;
        else
            l = ceil(h(j + 1) / ds);
        end

        % Beam update at each node
        if i == l
            w(i, j + 1) = S21 * wssss + C + 2 * w(i, j) - w(i, j - 1);
        else
            w(i, j + 1) = S2 * wssss + C + 2 * w(i, j) - w(i, j - 1);
        end

        % This is how w3 was computed in the original (overwritten each i)
        w3(j + 1) = w(i, j + 1) - w(1, j + 1);

        % Boundary conditions
        wsssl = (-2 * w(n, j) + 3 * w(n - 1, j) - w(n - 2, j));
        S3    = (EI * (dt^2)) / (mk * 2 * ds^3);

        w(2, j + 1)     = w(1, j + 1);
        w(n, j + 1)     = 2 * w(n, j) - w(n, j - 1) + S3 * wsssl;
        w(n - 1, j + 1) = (w(n, j + 1) + w(n - 2, j + 1)) / 2;

        %% Controller
        switch flag
            case 0
                F1(1:2000)  = 10;
                F1(2001:4000) = -10;
                F2(1:2000)  = 10;
                F2(2001:4000) = -10;
                check2 = zeros(1, r); 

            case 1
                % ROBUST-PD
                k1 = 10;
                kx = 20;
                kn = 1.5;

                dx1     = (w(1, j + 1) - w(1, j)) / dt;
                F1(j + 2) = -k1 * (w(1, j + 1) - sp1) - kx * dx1;

                S2_ = (w(l + 1, j + 1) - 2 * w(l, j + 1) + w(l - 1, j + 1)) / (ds^3);
                dx2 = ((w(l + 1, j + 1) - w(l - 1, j + 1)) / 2) * ds;
                dtx2 = (h(j + 1) - h(j)) / dt;

                F2(j + 2) = mh * g + EI * S2_ * dx2 - (h(j + 1) - sp2) - kn * dtx2;

            case 2
                % barrier
                k1   = 20;
                kmax = 0.004;
                k2   = 0.01;
                k3   = 30;
                k0   = 10;
                kc   = 0.001;
                kn   = 1.5;

                z   = w(l, j + 1);
                dwl = (w(1, j + 1) - w(1, j)) / dt;

                S1 = kc / ((kmax^2 - z^2) * 2.3) + k2;
                F1(j + 2) = -k1 * (w(1, j + 1) - sp1) - k3 * dwl - S1 * kmax * k0 * dwl - S1 * kmax;

                S2_  = (w(l + 1, j + 1) - 2 * w(l, j + 1) + w(l - 1, j + 1)) / (ds^3);
                dx2  = ((w(l + 1, j + 1) - w(l - 1, j + 1)) / 2) * ds;
                dtx2 = (h(j + 1) - h(j)) / dt;

                F2(j + 2) = mh * g + EI * S2_ * dx2 - (h(j + 1) - sp2) - kn * dtx2;

            case 3
                % PD
                ki1 = 1;
                kd1 = 10;

                e(j + 1)  = (sp1 - w(1, j + 1));
                de        = (e(j + 1) - e(j)) / dt;
                F1(j + 2) = ki1 * (sp1 - w(1, j + 1)) + kd1 * de;

            case 4
                % disturbance observer
                a1 = 1;
                k1 = 2.5;
                k2 = 25;
                kn = 1.5;

                dwl = (w(1, j + 1) - w(1, j)) / dt;

                dz1      = -a1 * (F1(j + 1) - EI * wsss0) - a1 * (z1(j + 1) + a1 * mw * dwl);
                z1(j + 1) = dz1 * dt + z1(j);

                d1e       = z1(j + 1) + a1 * mw * dwl;
                F1(j + 2) = -k1 * dwl - d1e - k2 * (w(1, j + 1) - sp1);

                S2_  = (w(l + 1, j + 1) - 2 * w(l, j + 1) + w(l - 1, j + 1)) / (ds^3);
                dx2  = (w(l + 1, j + 1) - w(l - 1, j + 1)) / ds;
                dtx2 = (h(j + 1) - h(j)) / dt;

                F2(j + 2) = mh * g + EI * S2_ * dx2 - (h(j + 1) - sp2) - kn * dtx2;

            case 5
                % (kept from original)
                um = 15;
                vm = 60;
                c1 = 10;
                c3 = 13;
                c2 = 40;
                k1 = 0.12;
                kn = 1.5;

                dw1  = (w(1, j + 1) - w(1, j)) / dt;
                wsss = (w(3, j + 1) - 3 * w(2, j + 1) + 2 * w(1, j + 1)) / ds(j + 1)^3; 

                e(j + 1) = (sp1 - w(1, j + 1));
                de       = (e(j + 1) - e(j)) / dt;

                a1(j + 1) = -c1 * dw1 + 1/2 * EI * wsss + 1/2 * k1 * e(j + 1);
                d11       = dw1;
                d21       = F1(j + 1) - a1(j + 1);

                a2(j + 1) = -2 * (1 / mw) * d11 - (1 / mw) * c2 * d21 + (a1(j + 1) - a1(j)) / dt;
                d31(j + 1) = dF1(j + 1) - a2(j + 1); 

                da2     = (a2(j + 1) - a2(j)) / dt;
                U1(j + 2) = (-1 / mw) * c3 * d31(j + 1) - d21 + da2; 

                dv1       = (1 / (cosh(v1(j + 1) / vm))^2) * U1(j + 2);
                v1(j + 2) = dt * dv1 + v1(j + 1);

                dF1(j + 2) = vm * tanh(v1(j + 2) / vm);
                du1        = (1 / (cosh(u1(j + 1) / vm))^2) * dF1(j + 2);

                u1(j + 2) = du1 * dt + u1(j + 1);
                F1(j + 2) = um * tanh(u1(j + 2) / um);

                S2_  = (w(l + 1, j + 1) - 2 * w(l, j + 1) + w(l - 1, j + 1)) / (ds^3);
                dx2  = ((w(l + 1, j + 1) - w(l - 1, j + 1)) / 2) * ds;
                dtx2 = (h(j + 1) - h(j)) / dt;

                F2(j + 2) = mh * g + EI * S2_ * dx2 - (h(j + 1) - sp2) - kn * dtx2;

            case 6
                % (kept from original)
                kmax = 0.01;
                g1   = 0.1;
                k1   = 1;
                c1   = 1;
                k11  = 1;
                k2   = 1;

                wsss0_ = (w(3, j) - 3 * w(2, j) + 2 * w(1, j)); 
                wsssl_ = (-2 * w(n, j) + 3 * w(n - 1, j) - w(n - 2, j)); 

                dw1 = (w(1, j + 1) - w(1, j)) / dt;
                z(j + 1) = w(n, j + 1) - w(1, j + 1); 

                a = log((2 * kmax^2) / (kmax^2 - z(j + 1)^2)) / log(2.71);
                dz = (z(j + 1) - z(j)) / dt;

                s1 = EI * wsss0_ / ds(j)^3; 

                v1(j + 2) = -k1 * dz - s1 - k11 * dz / a - mw * dz / a * z(j + 1) * dz / (kmax^2 - z(j + 1)^2) ...
                    - (k2 * (w(1, j) - sp1)) / a - mw * EI * wsssl_ / (ds(j)^3 * mk(j) * a);

                p1(j + 2)  = (g1 * dz * v1(j + 2) * a - p1(j + 1)) * dt + p1(j + 1);
                Fc1(j + 2) = p1(j + 2) * v1(j + 2);
                F1(j + 2)  = c1 * Fc1(j + 2);

            case 7
                % ADRC + ZVD xe con / xe nâng (kept from original)
                w1_eso(:, j + 2) = A_ESO1 * w1_eso(:, j + 1) + B_ESO1 * F1(j + 1) + Lc1 * w(1, j + 1);
                F1(j + 2) = (Kp1 * (sp1 * sp1_is(j) - w1_eso(1, j + 2)) - Kd1 * w1_eso(2, j + 2) - w1_eso(3, j + 2)) / b01;

                w2_eso(:, j + 2) = A_ESO2 * w2_eso(:, j + 1) + B_ESO2 * (F2(j + 1) - mh * g) + Lc2 * h(j + 1);
                F2(j + 2) = (Kp2 * (sp2 - w2_eso(1, j + 2)) - Kd2 * w2_eso(2, j + 2) - w2_eso(3, j + 2)) / b02 + mh * g;

            otherwise
                % No controller update
        end
    end
end

%% Plotting
x = 0 : dt : tmax;
y = 0 : ds : L;

figure(1)
subplot(2, 2, 1)
grid on; hold on;
plot(x, w(1, :), 'b', LineWidth=1);
xlabel('Thời gian(s)');
ylabel('Vị trí xe con(m)');
title('Vị trí xe con');
ylim([0 2]);

subplot(2, 2, 2)
hold on; grid on;

dolac = zeros(1, r);
for j = 1 : r
    dolac(j) = w(n, j) - w(1, j);
end

plot(x, dolac, 'b', LineWidth=1);
title('Độ lắc điểm cuối');
xlabel('Thời gian(s)');
ylabel('m');

subplot(2, 2, 3)
plot(x, h, 'r', LineWidth=1);
grid on;
title('Vị trí xe nâng');
xlabel('Thời gian(s)');
ylabel('m');

subplot(2, 2, 4)
grid on;
plot(x, w3, LineWidth=1);
title('Độ lắc xe nâng');
xlabel('Thời gian(s)');
ylabel('m');

figure(2)
subplot(1, 2, 1)
plot(x, F1(1:r), LineWidth=1.5);
title('F1');
grid on;

check = dolac - w3;
subplot(1, 2, 2)
plot(x, F2(1:r), LineWidth=1);
title('F2');
grid on;

figure(3)
% Phân tích phổ
Fs   = 1 / dt;
time = 0 : dt : (r - 3) * dt;
l1   = r;

freq = 0 : (1 / time(end)) : Fs/2 - (1 / time(end));

for i = 1 : n
    fft_w = fft(w(i, :) - w(1, :), l1) * (2 / l1);
    abs_w = abs(fft_w);

    subplot(ceil(n / 2), 2, i);
    plot(freq, abs_w(1:length(freq)), 'LineWidth', 0.8);
    hold on; grid on;

    set(gca, 'FontName', 'Times New Roman', 'FontSize', 9);
    xlabel('Frequency(Hz)');
    ylabel('Amplitude(m)');

    title(['FFT của thanh tại điểm ', num2str(i), ' khi xe con ở vị trí ', num2str(l)]);
end
