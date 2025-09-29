% stacker_crane_two_forces.m
% Mô phỏng stacker crane với 2 lực điều khiển:
%   F1 -> carriage (xe con) theo phương ngang
%   F2 -> elevator / moving mass (di chuyển dọc dầm)
% Dầm: Euler-Bernoulli (5-point FD), tích phân thời gian Newmark-beta

clear; close all; clc;

%% ---------- Tham số vật lý & số ----------
L = 1.0;            % chiều dài dầm (m)
E = 1e4;            % Young's modulus
I = 1e-4;           % second moment of area
rho = 1.0;          % density (mass per unit volume)
A = 1.0;            % diện tích tiết diện

% Masses
M_c = 50;           % carriage mass (kg)
c_c = 50;           % carriage viscous damping (N s/m)

M_l = 7.5;          % elevator / moving load mass (kg)
c_l = 20;           % elevator damping (N s/m)

% beam discretization
N = 50;
dx = L / N;

%% ---------- Beam matrices ----------
K = assemble_K(N, E, I, dx);
Mb = rho * A * dx * eye(N);   % lumped mass of beam nodes

%% ---------- Time & Newmark ----------
beta = 1/4; gamma = 1/2;
dt = 0.005;
T_end = 10.0;
t_arr = (0:dt:T_end)';
nt = numel(t_arr);
F1 = zeros(nt,1);
F2 = zeros(nt,1);
%% ---------- Control forces (hằng số) ----------
% F1_const = 10;      % N, force on carriage
% F2_const = 10;       % N, force on elevator (up/down along beam axis)
% 
% F1 = F1_const * ones(nt,1);
% F2 = F2_const * ones(nt,1);
%% ---------- Control gains & setpoints ----------
Kp1 = 500;  % N/m
Kd1 = 50;   % N s/m
x_c_ref = 0.8*L;  % target position for carriage

Kp2 = 200;  % N/m
Kd2 = 20;   % N s/m
s_ref = 0.75*L;   % target position for moving mass

%% ---------- State arrays ----------
% Carriage
xc = zeros(nt,1);
xdot = zeros(nt,1);
xddot = zeros(nt,1);

% Elevator (moving load) states: s = position along beam (m)
s = zeros(nt,1);        % initial position s(1)=0 -> you can set initial
sdot = zeros(nt,1);
sddot = zeros(nt,1);
% Initialize elevator starting position (ví dụ bắt đầu ở 0.25L)
s(1) = 0.25 * L;

% Beam states
w = zeros(N,1);
wdot = zeros(N,1);
wddot = zeros(N,1);
W_rel = zeros(nt, N+1); % col 1 = x=0 -> 0

%% ---------- Precompute Newmark constants ----------
a1 = 1/(beta*dt^2);
a2 = 1/(beta*dt);
a3 = 1/(2*beta)-1;

%% ---------- Time loop ----------
for k = 1:nt
   F1(k) = Kp1*(x_c_ref - xc(k)) - Kd1*xdot(k);
    F2(k) = Kp2*(s_ref - s(k)) - Kd2*sdot(k);
    %% --- 1) Carriage integrator (explicit Euler) ---
    % dynamic: M_c * xddot = F1(k) - c_c * xdot
    if k == 1
        xddot(k) = (F1(k) - c_c * xdot(k)) / M_c;
        xdot(k) = xdot(k) + xddot(k) * dt;
        xc(k) = xc(k) + xdot(k) * dt;
    else
        xddot(k) = (F1(k) - c_c * xdot(k-1)) / M_c;
        xdot(k) = xdot(k-1) + xddot(k) * dt;
        xc(k) = xc(k-1) + xdot(k) * dt;
    end

    %% --- 2) Elevator integrator (explicit Euler) ---
    % dynamic: M_l * sddot = F2(k) - c_l * sdot(k)  (we ignore direct beam reaction for simplicity)
    if k == 1
        sddot(k) = (F2(k) - c_l * sdot(k)) / M_l;
        sdot(k) = sdot(k) + sddot(k) * dt;
        s(k) = s(k) + sdot(k) * dt;
    else
        sddot(k) = (F2(k) - c_l * sdot(k-1)) / M_l;
        sdot(k) = sdot(k-1) + sddot(k) * dt;
        s(k) = s(k-1) + sdot(k) * dt;
    end

    % Ensure s stays inside beam [0, L]
    if s(k) < 0, s(k) = 0; sdot(k) = 0; end
    if s(k) > L, s(k) = L; sdot(k) = 0; end

    %% --- 3) Beam: base acceleration is carriage acceleration xddot(k) ---
    a_base_k = xddot(k);

    % Build moving mass matrix Mm at time k using current s(k)
    Mm = zeros(N,N);
    idx = floor(s(k) / dx);            % 0-based
    alpha = (s(k) - idx*dx) / dx;
    if idx >= 0 && idx < N
        Mm(idx+1, idx+1) = Mm(idx+1, idx+1) + M_l * (1 - alpha);
    end
    if (idx+1) >= 0 && (idx+1) < N
        Mm(idx+2, idx+2) = Mm(idx+2, idx+2) + M_l * alpha;
    end
    M = Mb + Mm;

    % Save relative displacement (left end x=0 treated 0)
    W_rel(k,1) = 0;
    W_rel(k,2:end) = w(:)';

    if k == nt
        break;
    end

    %% --- 4) Predictor for next step carriage acceleration & elevator accel (explicit) ---
    % Predictor uses current velocities + forces at k+1 (here F1(k+1), F2(k+1))
    % Carriage predict
    xddot_next = (F1(min(k+1,nt)) - c_c * xdot(k)) / M_c;
    % Elevator predict
    sddot_next = (F2(min(k+1,nt)) - c_l * sdot(k)) / M_l;

    a_base_next = xddot_next;

    % Predict s_next as well (for mass distribution at k+1) - simple Euler predictor
    s_pred = s(k) + sdot(k) * dt + 0.5 * sddot_next * dt^2;
    % bound to [0,L]
    if s_pred < 0, s_pred = 0; end
    if s_pred > L, s_pred = L; end

    %% --- 5) Build M_next using s_pred (approx s at k+1) ---
    Mm_next = zeros(N,N);
    idx_n = floor(s_pred / dx);
    alpha_n = (s_pred - idx_n*dx) / dx;
    if idx_n >= 0 && idx_n < N
        Mm_next(idx_n+1, idx_n+1) = Mm_next(idx_n+1, idx_n+1) + M_l * (1 - alpha_n);
    end
    if (idx_n+1) >= 0 && (idx_n+1) < N
        Mm_next(idx_n+2, idx_n+2) = Mm_next(idx_n+2, idx_n+2) + M_l * alpha_n;
    end
    M_next = Mb + Mm_next;

    %% --- 6) External force on beam due to base acceleration at n+1 ---
    f_np1 = - M_next * (a_base_next * ones(N,1));

    %% --- 7) Newmark implicit solve for beam ---
    K_eff = K + a1 * M_next;
    RHS = f_np1 + M_next * ( a1*w + a2*wdot + a3*wddot );

    w_new = K_eff \ RHS;

    % Update beam accel and vel
    wddot_new = a1*(w_new - w) - a2*wdot - a3*wddot;
    wdot_new = wdot + dt*((1-gamma)*wddot + gamma*wddot_new);

    % Assign
    w = w_new; wdot = wdot_new; wddot = wddot_new;

    %% --- 8) finalize predicted carriage & elevator states for next step ---
    xdot(k+1) = xdot(k) + xddot_next * dt;
    xc(k+1) = xc(k) + xdot(k+1) * dt;
    xddot(k+1) = xddot_next;

    sdot(k+1) = sdot(k) + sddot_next * dt;
    s(k+1) = s(k) + sdot(k+1) * dt;
    sddot(k+1) = sddot_next;

    % enforce bounds on s
    if s(k+1) < 0, s(k+1) = 0; sdot(k+1) = 0; end
    if s(k+1) > L, s(k+1) = L; sdot(k+1) = 0; end
end

%% ---------- Simple checks ----------
assert(~any(isnan(W_rel(:))), 'NaN in W_rel');
assert(~any(isinf(W_rel(:))), 'Inf in W_rel');

%% ---------- Plot results ----------
figure('Units','normalized','Position',[0.05 0.05 0.85 0.75]);

subplot(3,1,1);
plot(t_arr, xc, 'LineWidth',1.4);
xlabel('t (s)'); ylabel('x_{carriage} (m)');
title('Vị trí xe con theo thời gian (carriage position)');
grid on;

subplot(3,1,2);
plot(t_arr, s, 'LineWidth',1.2);
xlabel('t (s)'); ylabel('s(t) (m)');
title('Vị trí xe nâng (moving mass) theo thời gian');
ylim([0 L]); grid on;

subplot(3,1,3);
tip_rel = W_rel(:, end);
plot(t_arr, tip_rel, 'LineWidth',1.2);
xlabel('t (s)'); ylabel('Tip deflection (m)');
title('Dao động đỉnh dầm (tip deflection)');
grid on;

sgtitle('Stacker Crane Simulation (F1 & F2 constant)');

%% ========== Local functions ==========
function K = assemble_K(dof, E, I, dx)
    K = zeros(dof,dof);
    offsets = [-2 -1 0 1 2]; coeffs = [1 -4 6 -4 1];
    for j = 1:dof
        i = j;
        for s = 1:length(offsets)
            p = i + offsets(s);
            c = coeffs(s);
            if p == 0
                continue;
            elseif p < 0
                K(j,1) = K(j,1) + c;
            elseif p >= 1 && p <= dof
                K(j,p) = K(j,p) + c;
            elseif p == dof + 1
                K(j,dof) = K(j,dof) + c*2;
                if dof-1 >= 1, K(j,dof-1) = K(j,dof-1) - c; end
            elseif p == dof + 2
                K(j,dof) = K(j,dof) + c*4;
                if dof-1 >= 1, K(j,dof-1) = K(j,dof-1) - c*4; end
                if dof-2 >= 1, K(j,dof-2) = K(j,dof-2) + c; end
            end
        end
    end
    K = (E * I / dx^4) * K;
end
