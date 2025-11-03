% stacker_crane_sim_v2.m
% Fixed version: removes "value assigned might be unused" warnings
% Save as stacker_crane_sim_v2.m and run stacker_crane_sim_v2

function stacker_crane_sim_v2
    close all; clc;

    %% ---------- System parameters ----------
    p.M  = 1.5;
    p.m1 = 0.075;
    p.m2 = 0.15;
    p.m3 = 0.03;
    p.L  = 0.4;
    p.k  = 3.26;
    p.g  = 9.81;

    % friction/disturbance demo values
    p.frx = 0; p.epsx = 1e-3; p.krx = 0;
    p.fry = 0; p.epsy = 1e-3; p.kry = 0;
    p.D1 = 0; p.D2 = 0;

    %% ---------- Simulation settings ----------
    Tsim = 6.0;
    tspan = [0 Tsim];

    % initial state: [x; xdot; y; ydot; th1; dth1; th2; dth2; th3; dth3; th4; dth4]
    state0 = [0; 0; 0.05; 0;  0.01; 0; 0.01; 0; 0.01; 0; 0.01; 0];

    %% ---------- Reference / Gains ----------
    x0_ref = 0; xd = 0.5;
    y0_ref = 0.05; yd = 0.35;
    Tmove = 2.0;

    Kpx = 200; Kdx = 30;
    Kpy = 150; Kdy = 20;

    forceFun.Fa1 = @(t, state) ff_force_trolley(t, state, p, x0_ref, xd, Tmove, Kpx, Kdx);
    forceFun.Fa2 = @(t, state) ff_force_hoist (t, state, p, y0_ref, yd, Tmove, Kpy, Kdy);

    %% ---------- Integrate ----------
    opts = odeset('RelTol',1e-6, 'AbsTol',1e-8);
    try
        [T, S] = ode45(@(t, st) crane_ode(t, st, p, forceFun), tspan, state0, opts);
    catch ME
        fprintf('Error running ode45:\n%s\n', ME.message);
        rethrow(ME);
    end

    %% ---------- Build reference vectors (vectorized, avoids temp loop vars) ----------
    % use arrayfun so no loop variable remains unused and triggers warnings
    xref_vec = arrayfun(@(tt) traj_cos_profile(x0_ref, xd, Tmove, tt), T(:));
    yref_vec = arrayfun(@(tt) traj_cos_profile(y0_ref, yd, Tmove, tt), T(:));

    %% ---------- Plot ----------
    figure('Name','Stacker crane simulation (v2)');
    subplot(3,1,1);
    plot(T, S(:,1),'-b','LineWidth',1.3); hold on;
    plot(T, S(:,3),'-r','LineWidth',1.1);
    plot(T, xref_vec,'--k','LineWidth',1.0);
    plot(T, yref_vec,'--m','LineWidth',1.0);
    legend('x (trolley)','y (hoist)','x_{ref}','y_{ref}','Location','best');
    xlabel('t (s)'); ylabel('m'); grid on; title('Trolley & Hoist');

    subplot(3,1,2);
    plot(T, S(:,5),'LineWidth',1.0); hold on;
    plot(T, S(:,7),'LineWidth',1.0);
    plot(T, S(:,9),'LineWidth',1.0);
    plot(T, S(:,11),'LineWidth',1.0);
    legend('\theta_1','\theta_2','\theta_3','\theta_4','Location','best');
    xlabel('t (s)'); ylabel('rad'); grid on; title('Modal angles');

    subplot(3,1,3);
    plot(T, S(:,2),'-','LineWidth',1.0); hold on;
    plot(T, S(:,4),'-','LineWidth',1.0);
    legend('xdot','ydot','Location','best');
    xlabel('t (s)'); ylabel('m/s'); grid on; title('Velocities');
end

%% ---------------- ODE right-hand side ----------------
function dstate = crane_ode(t, state, p, forceFun)
    if numel(state) ~= 12
        error('State vector must have 12 elements.');
    end

    % Unpack
    X   = state(1);  Xd = state(2);
    Y   = state(3);  Yd = state(4);
    th1 = state(5);  dth1 = state(6);
    th2 = state(7);  dth2 = state(8);
    th3 = state(9);  dth3 = state(10);
    th4 = state(11); dth4 = state(12);

    % Forces
    Fa1 = forceFun.Fa1(t, state);
    Fa2 = forceFun.Fa2(t, state);

    % Friction
    Fr1 = p.frx * tanh(Xd / max(p.epsx,1e-6)) - p.krx * abs(Xd) * Xd;
    Fr2 = p.fry * tanh(Yd / max(p.epsy,1e-6)) - p.kry * abs(Yd) * Yd;

    % Disturbances
    Fd1 = p.D1; Fd2 = p.D2;

    F1 = Fa1 - Fr1 - Fd1;
    F2 = Fa2 - Fr2 - Fd2;

    % Build system and solve
    [A, b] = build_Ab(state, p, F1, F2);

    if rcond(A) < 1e-12
        A = A + 1e-9 * eye(6);
    end

    acc = A \ b;

    dstate = zeros(12,1);
    dstate(1) = Xd;    dstate(2) = acc(1);
    dstate(3) = Yd;    dstate(4) = acc(2);
    dstate(5) = dth1;  dstate(6) = acc(3);
    dstate(7) = dth2;  dstate(8) = acc(4);
    dstate(9) = dth3;  dstate(10)= acc(5);
    dstate(11)= dth4;  dstate(12)= acc(6);
end

%% --------------- Build A and b ----------------
function [A, b] = build_Ab(state, p, F1, F2)
    % Unpack
    Y   = state(3);  Yd = state(4);
    th1 = state(5);  dth1 = state(6);
    th2 = state(7);  dth2 = state(8);
    th3 = state(9);  dth3 = state(10);
    th4 = state(11); dth4 = state(12);

    M = p.M; m1 = p.m1; m2 = p.m2; m3 = p.m3; L = p.L; k = p.k; g = p.g;

    % clamp Y
    Y = min(max(Y,0), L);

    % trig terms
    c1 = cos(th1); s1 = sin(th1);
    c2 = cos(th2); s2 = sin(th2);
    c3 = cos(th3); s3 = sin(th3);
    c4 = cos(th4); s4 = sin(th4);


    c12 = cos(th1 - th2); s12 = sin(th1 - th2);
    c13 = cos(th1 - th3); s13 = sin(th1 - th3);
    c14 = cos(th1 - th4); s14 = sin(th1 - th4);
    c23 = cos(th2 - th3); s23 = sin(th2 - th3);
    c24 = cos(th2 - th4); s24 = sin(th2 - th4);
    c34 = cos(th3 - th4); s34 = sin(th3 - th4);

    % region offset
    if (Y >= 0) && (Y < L/4)
        Ybar = Y*s1;
    elseif (Y >= L/4) && (Y < L/2)
        Ybar = (Y - L/4)*s2+L/4*s1;
    elseif (Y >= L/2) && (Y < 3*L/4)
        Ybar = (Y - L/2)*s3+L/4*(s1+s2);
    else
        Ybar = (Y - 3*L/4)*s4+L/4*(s1+s2+s3);
    end

    coef1_th1 = L * (7*m1/32 + m3/4) + m2 * Ybar;% coefficient=hệ số 
    coef1_th2 = (5*m1/32 + m3/4) * L;
    coef1_th3 = (3*m1/32 + m3/4) * L;
    coef1_th4 = (m1/32 + m3/4) * L;

    A = zeros(6,6); b = zeros(6,1);

    % Phương trình  (1)
    A(1,1) = M + m1 + m2 + m3;
    A(1,2) = m2 * s1;
    A(1,3) = coef1_th1 * c1;
    A(1,4) = coef1_th2 * c2;
    A(1,5) = coef1_th3 * c3;
    A(1,6) = coef1_th4 * c4;
    b(1) = F1 + coef1_th1 * s1 * dth1^2 + coef1_th2 * s2 * dth2^2 ...
           + coef1_th3 * s3 * dth3^2 + coef1_th4 * s4 * dth4^2 ...
           - 2 * m2 * c1 * Yd * dth1;

    % Phương trình (2)
    A(2,1) = m2 * s1;
    A(2,2) = m2;
    b(2) = F2 + m2 * Ybar * dth1^2 + m2 * g * c1;

    % Phương trình  (3)
    A(3,1) = (m2 * Ybar + L * (7*m1/32 + m3/4)) * c1;
    A(3,3) = m2 * Ybar^2 + L^2 * (13*m1/256 + m3/16);
    A(3,4) = (5*m1/128 + m3/16) * L^2 * c12;
    A(3,5) = (3*m1/128 + m3/16) * L^2 * c13;
    A(3,6) = (m1/128 + m3/16) * L^2 * c14;
    NL3 = (5*m1/128 + m3/16) * L^2 * s12 * dth2^2 ...
        + (3*m1/128 + m3/16) * L^2 * s13 * dth3^2 ...
        + (m1/128 + m3/16) * L^2 * s14 * dth4^2 ...
        + 2 * m2 * Ybar * Yd * dth1 + k * th1 ...
        - m2 * g * s1 * Ybar - (7*m1/32 + m3/4) * g * L * s1;
    b(3) = - NL3;

    % Phương trình  (4)
    A(4,1) = (5*m1/32 + m3/4) * L * c2;
    A(4,3) = (5*m1/128 + m3/16) * L^2 * c12;
    A(4,4) = (9*m1/256 + m3/16) * L^2;
    A(4,5) = (3*m1/128 + m3/16) * L^2 * c23;
    A(4,6) = (m1/128 + m3/16) * L^2 * c24;
    NL4 = - (5*m1/128 + m3/16) * L^2 * s12 * dth1^2 ...
        + (3*m1/128 + m3/16) * L^2 * s23 * dth3^2 ...
        + (m1/128 + m3/16) * L^2 * s24 * dth4^2 ...
        + k * th2 - (5*m1/32 + m3/4) * g * L * s2;
    b(4) = - NL4;

    % Phương trình  (5)
    A(5,1) = (3*m1/32 + m3/4) * L * c3;
    A(5,3) = (3*m1/128 + m3/16) * L^2 * c13;
    A(5,4) = (3*m1/128 + m3/16) * L^2 * c23;
    A(5,5) = (5*m1/256 + m3/16) * L^2;
    A(5,6) = (m1/128 + m3/16) * L^2 * c34;
    NL5 = - (3*m1/128 + m3/16) * L^2 * s13 * dth1^2 ...
        - (3*m1/128 + m3/16) * L^2 * s23 * dth2^2 ...
        + (m1/128 + m3/16) * L^2 * s34 * dth4^2 ...
        + k * th3 - (3*m1/32 + m3/4) * g * L * s3;
    b(5) = - NL5;

    % Phương trình  (6)
    A(6,1) = (m1/32 + m3/4) * L * c4;
    A(6,3) = (m1/128 + m3/16) * L^2 * c14;
    A(6,4) = (m1/128 + m3/16) * L^2 * c24;
    A(6,5) = (m1/128 + m3/16) * L^2 * c34;
    A(6,6) = (m1/256 + m3/16) * L^2;
    NL6 = - (m1/128 + m3/16) * L^2 * s14 * dth1^2 ...
        - (m1/128 + m3/16) * L^2 * s24 * dth2^2 ...
        - (m1/128 + m3/16) * L^2 * s34 * dth3^2 ...
        + k * th4 - (m1/32 + m3/4) * g * L * s4;
    b(6) = - NL6;
end

%% ---------- helpers (FF+PD & traj) ----------
function F = ff_force_trolley(t, state, p, x0, xf, Tmove, Kp, Kd)
    [pos, vel, acc] = traj_cos_profile_deriv(x0, xf, Tmove, t);
    curr_x = state(1); curr_v = state(2);
    Mtot = p.M + p.m1 + p.m2 + p.m3;
    F = Mtot * acc + Kp * (pos - curr_x) + Kd * (vel - curr_v);
end

function F = ff_force_hoist(t, state, p, y0, yf, Tmove, Kp, Kd)
    [pos, vel, acc] = traj_cos_profile_deriv(y0, yf, Tmove, t);
    curr_y = state(3); curr_v = state(4);
    Mhoist = p.m2; if Mhoist < 1e-6, Mhoist = 0.01; end
    F = Mhoist * acc + Kp * (pos - curr_y) + Kd * (vel - curr_v);
end

function s = traj_cos_profile(x0, xf, T, t)
    if t <= 0, s = x0;
    elseif t >= T, s = xf;
    else A=(xf-x0)/2; B=pi/T; s = x0 + A*(1-cos(B*t));
    end
end

function [pos, vel, acc] = traj_cos_profile_deriv(x0, xf, T, t)
    A = (xf - x0) / 2; B = pi / T;
    if t <= 0, pos = x0; vel = 0; acc = 0;
    elseif t >= T, pos = xf; vel = 0; acc = 0;
    else pos = x0 + A*(1 - cos(B*t)); vel = A * B * sin(B*t); acc = A * B^2 * cos(B*t);
    end
end
