%% Stacker Crane Dynamics Simulation with Separate Pendulum Models
% Four separate ODE functions for each payload position case

clear; close all; clc;

%% System Parameters
Parameters = struct();
Parameters.L = 0.4;      % Total mast length (m)
Parameters.M = 1.5;      % Trolley mass (kg)
Parameters.m1 = 0.075;   % First pendulum mass (kg)
Parameters.m2 = 0.15;    % Payload mass (kg)
Parameters.m3 = 0.03;    % Tip mass (kg)
Parameters.g = 9.81;     % Gravitational acceleration (m/s^2)
Parameters.k = 3.26;     % Virtual spring stiffness (N*m/rad)

L = Parameters.L;

%% Initial Conditions
% State vector: [x, y, th1, th2, th3, th4, xd, yd, th1d, th2d, th3d, th4d]
state0 = [0;           % x: Initial trolley position (m)
          0.05;        % y: Initial payload height (m)
          0.01;        % th1: Initial angle of pendulum 1 (rad)
          0.01;        % th2: Initial angle of pendulum 2 (rad)
          0.01;        % th3: Initial angle of pendulum 3 (rad)
          0.01;        % th4: Initial angle of pendulum 4 (rad)
          0;           % xd: Initial trolley velocity (m/s)
          0;           % yd: Initial payload velocity (m/s)
          0;           % th1d: Initial angular velocity 1 (rad/s)
          0;           % th2d: Initial angular velocity 2 (rad/s)
          0;           % th3d: Initial angular velocity 3 (rad/s)
          0];          % th4d: Initial angular velocity 4 (rad/s)

%% Simulation Parameters
t_end = 10;       % Total simulation time (s)
dt = 0.01;        % Time step (s)

%% Control Forces
t_span = 0:dt:t_end;
n_steps = length(t_span);

% Preallocate force arrays
F1_array = zeros(n_steps, 1);  % Horizontal force on trolley
F2_array = zeros(n_steps, 1);  % Vertical force on payload

% Define trapezoidal velocity profile for trolley
t_accel = 2;
t_decel = 8;
F_max = 2;

for i = 1:n_steps
    t = t_span(i);
    if t < t_accel
        F1_array(i) = F_max;
    elseif t >= t_accel && t < t_decel
        F1_array(i) = 0;
    else
        F1_array(i) = -F_max;
    end
    
    % Payload hoisting command
    if t < 6
        F2_array(i) = 0.3;
    else
        F2_array(i) = 0;
    end
end

%% Initialize Storage Arrays
states = zeros(n_steps, 12);
states(1, :) = state0';
w_tip = zeros(n_steps, 1);

%% Main Simulation Loop with Event-Based Model Switching
fprintf('Starting simulation with event detection...\n');

% Track which model we're using
t_current = 0;
state_current = state0;
y_current = state_current(2);

% Determine initial model
if y_current >= 0 && y_current < L/4
    current_model = 1;
elseif y_current >= L/4 && y_current < L/2
    current_model = 2;
elseif y_current >= L/2 && y_current < 3*L/4
    current_model = 3;
else
    current_model = 4;
end

% Storage for results
t_all = [];
states_all = [];
segment = 1;

while t_current < t_end
    fprintf('  Segment %d: Using Model %d at t=%.2fs, y=%.3fm\n', ...
            segment, current_model, t_current, y_current);
    
    % Set up event detection for boundary crossings
    options = odeset('Events', @(t, state) boundary_event(t, state, L), ...
                     'RelTol', 1e-6, 'AbsTol', 1e-8);
    
    % Interpolate forces for current segment
    Forces_func = @(t) [interp1(t_span, F1_array, t, 'linear', 'extrap'); ...
                        interp1(t_span, F2_array, t, 'linear', 'extrap')];
    
    % Select appropriate model
    switch current_model
        case 1
            ode_func = @(t, state) first_pendulum_event(t, state, Parameters, Forces_func);
        case 2
            ode_func = @(t, state) second_pendulum_event(t, state, Parameters, Forces_func);
        case 3
            ode_func = @(t, state) third_pendulum_event(t, state, Parameters, Forces_func);
        case 4
            ode_func = @(t, state) fourth_pendulum_event(t, state, Parameters, Forces_func);
    end
    
    % Solve until next boundary or end time
    [t_seg, state_seg, te, ye, ie] = ode45(ode_func, [t_current, t_end], state_current, options);
    
    % Store results
    t_all = [t_all; t_seg];
    states_all = [states_all; state_seg];
    
    % Check if we hit a boundary
    if ~isempty(ie)
        % Event detected - update state and determine new model
        t_current = te(end);
        state_current = ye(end, :)';
        y_current = state_current(2);
        
        % Determine which model to use next
        if y_current >= 0 && y_current < L/4
            current_model = 1;
        elseif y_current >= L/4 && y_current < L/2
            current_model = 2;
        elseif y_current >= L/2 && y_current < 3*L/4
            current_model = 3;
        else
            current_model = 4;
        end
        
        segment = segment + 1;
        
        % Safety check
        if segment > 100
            warning('Too many model switches, stopping');
            break;
        end
    else
        % Reached end time without hitting boundary
        fprintf('  Simulation completed at t=%.2fs\n', t_end);
        break;
    end
end

fprintf('Total segments: %d\n', segment);

% Interpolate to uniform time grid for plotting
states = interp1(t_all, states_all, t_span);

% Calculate tip displacement
for i = 1:length(t_span)
    y_val = states(i, 2);
    th1 = states(i, 3);
    th2 = states(i, 4);
    th3 = states(i, 5);
    th4 = states(i, 6);
    
    if y_val >= 0 && y_val < L/4
        w_tip(i) = y_val * sin(th1);
    elseif y_val >= L/4 && y_val < L/2
        w_tip(i) = L/4*sin(th1) + (y_val - L/4)*sin(th2);
    elseif y_val >= L/2 && y_val < 3*L/4
        w_tip(i) = L/4*(sin(th1) + sin(th2)) + (y_val - L/2)*sin(th3);
    else
        w_tip(i) = L/4*(sin(th1) + sin(th2) + sin(th3)) + (y_val - 3*L/4)*sin(th4);
    end
end

%% Plotting Results
figure('Position', [100, 100, 1200, 800]);

% Plot 1: Pendulum Angles
subplot(2, 2, 1);
plot(t_span, rad2deg(states(:, 3)), 'LineWidth', 1.5); hold on;
plot(t_span, rad2deg(states(:, 4)), 'LineWidth', 1.5);
plot(t_span, rad2deg(states(:, 5)), 'LineWidth', 1.5);
plot(t_span, rad2deg(states(:, 6)), 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Angle (degrees)');
title('Pendulum Angles');
legend('\theta_1', '\theta_2', '\theta_3', '\theta_4', 'Location', 'best');

% Plot 2: Tip Displacement
subplot(2, 2, 2);
plot(t_span, w_tip*100, 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
grid on;
xlabel('Time (s)');
ylabel('Tip Displacement (cm)');
title('Mast Tip Lateral Displacement');

% Plot 3: Trolley Position
subplot(2, 2, 3);
plot(t_span, states(:, 1), 'LineWidth', 1.5, 'Color', [0.47 0.67 0.19]);
grid on;
xlabel('Time (s)');
ylabel('Position (m)');
title('Trolley Horizontal Position');

% Plot 4: Payload Vertical Position
subplot(2, 2, 4);
plot(t_span, states(:, 2)*100, 'LineWidth', 1.5, 'Color', [0.49 0.18 0.56]);
hold on;
yline(L*100/4, '--r', 'L/4', 'LineWidth', 1.5);
yline(L*100/2, '--r', 'L/2', 'LineWidth', 1.5);
yline(3*L*100/4, '--r', '3L/4', 'LineWidth', 1.5);
yline(L*100, '--r', 'L', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Height (cm)');
title('Payload Vertical Position');
legend('Payload', 'Boundaries', 'Location', 'best');

sgtitle('Stacker Crane Dynamics Simulation', 'FontSize', 14, 'FontWeight', 'bold');

%% ========================================================================
%  SEPARATE DYNAMIC EQUATIONS FOR EACH PENDULUM POSITION
%  ========================================================================

%% Event Function to Detect Boundary Crossings
function [value, isterminal, direction] = boundary_event(~, state, L)
    % Extract payload position
    y = state(2);
    
    % Define the three boundaries where model switches occur
    boundaries = [L/4, L/2, 3*L/4];
    
    % Create events for each boundary (both upward and downward crossings)
    % We need 6 events total: 3 boundaries × 2 directions
    value = [y - L/4;      % Event 1: Crosses L/4 going up
             L/4 - y;      % Event 2: Crosses L/4 going down
             y - L/2;      % Event 3: Crosses L/2 going up
             L/2 - y;      % Event 4: Crosses L/2 going down
             y - 3*L/4;    % Event 5: Crosses 3L/4 going up
             3*L/4 - y];   % Event 6: Crosses 3L/4 going down
    
    % Stop integration when ANY boundary is crossed
    isterminal = ones(6, 1);
    
    % Direction: positive value means event triggers on zero crossing
    % All set to 1 means trigger when value crosses zero going positive
    direction = ones(6, 1);
end

%% Modified Pendulum Functions with Time-Varying Forces

%% First Pendulum - Event Version
function dstate = first_pendulum_event(t, state, P, Forces_func)
    % Get forces at current time
    Forces = Forces_func(t);
    
    % Call the original function
    dstate = first_pendulum(t, state, P, Forces);
end

%% Second Pendulum - Event Version
function dstate = second_pendulum_event(t, state, P, Forces_func)
    Forces = Forces_func(t);
    dstate = second_pendulum(t, state, P, Forces);
end

%% Third Pendulum - Event Version
function dstate = third_pendulum_event(t, state, P, Forces_func)
    Forces = Forces_func(t);
    dstate = third_pendulum(t, state, P, Forces);
end

%% Fourth Pendulum - Event Version
function dstate = fourth_pendulum_event(t, state, P, Forces_func)
    Forces = Forces_func(t);
    dstate = fourth_pendulum(t, state, P, Forces);
end

%% Helper Function: Compute All Trigonometric Terms
function trig = compute_trig(th1, th2, th3, th4)
    % Single, sine and cosine for each angle
    trig.s1 = sin(th1); trig.c1 = cos(th1);
    trig.s2 = sin(th2); trig.c2 = cos(th2);
    trig.s3 = sin(th3); trig.c3 = cos(th3);
    trig.s4 = sin(th4); trig.c4 = cos(th4);
    
    % Difference angles
    trig.s12 = sin(th1-th2); trig.c12 = cos(th1-th2);
    trig.s13 = sin(th1-th3); trig.c13 = cos(th1-th3);
    trig.s14 = sin(th1-th4); trig.c14 = cos(th1-th4);
    trig.s23 = sin(th2-th3); trig.c23 = cos(th2-th3);
    trig.s24 = sin(th2-th4); trig.c24 = cos(th2-th4);
    trig.s34 = sin(th3-th4); trig.c34 = cos(th3-th4);
end

%% First Pendulum Dynamics (0 <= y < L/4)
function dstate = first_pendulum(~, state, P, Forces)
    % Extract parameters
    L = P.L; M = P.M; m1 = P.m1; m2 = P.m2; m3 = P.m3; g = P.g; k = P.k;
    
    % Extract state variables
    y = state(2);
    th1 = state(3); th2 = state(4); th3 = state(5); th4 = state(6);
    yd = state(8);
    th1d = state(9); th2d = state(10); th3d = state(11); th4d = state(12);
    
    % Compute all trigonometric terms at once
    t = compute_trig(th1, th2, th3, th4);
    s1 = t.s1; c1 = t.c1; s2 = t.s2; c2 = t.c2;
    s3 = t.s3; c3 = t.c3; s4 = t.s4; c4 = t.c4;
    s12 = t.s12; c12 = t.c12; s13 = t.s13; c13 = t.c13;
    s14 = t.s14; c14 = t.c14; s23 = t.s23; c23 = t.c23;
    s24 = t.s24; c24 = t.c24; s34 = t.s34; c34 = t.c34;
    
    % Forces
    F1 = Forces(1);
    F2 = Forces(2);
    
    % Build mass matrix
    Mmat = zeros(6,6);
    Mmat(1,1) = M + m1 + m2 + m3;
    Mmat(1,2) = m2*s1;
    Mmat(1,3) = (7*m1/32*L + m3/4*L + m2*y)*c1;
    Mmat(1,4) = (5*m1/32 + m3/4)*L*c2;
    Mmat(1,5) = (3*m1/32 + m3/4)*L*c3;
    Mmat(1,6) = (m1/32 + m3/4)*L*c4;
    
    Mmat(2,1) = m2*s1;
    Mmat(2,2) = m2;
    
    Mmat(3,1) = (m2*y + 7*m1/32*L + m3/4*L)*c1;
    Mmat(3,3) = m2*y^2 + 13*m1/256*L^2 + m3/16*L^2;
    Mmat(3,4) = (5*m1/128 + m3/16)*L^2*c12;
    Mmat(3,5) = (3*m1/128 + m3/16)*L^2*c13;
    Mmat(3,6) = (m1/128 + m3/16)*L^2*c14;
    
    Mmat(4,1) = (5*m1/32 + m3/4)*L*c2;
    Mmat(4,3) = (5*m1/128 + m3/16)*L^2*c12;
    Mmat(4,4) = (9*m1/256 + m3/16)*L^2;
    Mmat(4,5) = (3*m1/128 + m3/16)*L^2*c23;
    Mmat(4,6) = (m1/128 + m3/16)*L^2*c24;
    
    Mmat(5,1) = (3*m1/32 + m3/4)*L*c3;
    Mmat(5,3) = (3*m1/128 + m3/16)*L^2*c13;
    Mmat(5,4) = (3*m1/128 + m3/16)*L^2*c23;
    Mmat(5,5) = (5*m1/256 + m3/16)*L^2;
    Mmat(5,6) = (m1/128 + m3/16)*L^2*c34;
    
    Mmat(6,1) = (m1/32 + m3/4)*L*c4;
    Mmat(6,3) = (m1/128 + m3/16)*L^2*c14;
    Mmat(6,4) = (m1/128 + m3/16)*L^2*c24;
    Mmat(6,5) = (m1/128 + m3/16)*L^2*c34;
    Mmat(6,6) = (m1/256 + m3/16)*L^2;
    
    % Build forcing vector
    Phi = zeros(6,1);
    Phi(1) = F1 + (7*m1/32*L + m3/4*L + m2*y)*s1*th1d^2 ...
             + (5*m1/32 + m3/4)*L*s2*th2d^2 ...
             + (3*m1/32 + m3/4)*L*s3*th3d^2 ...
             + (m1/32 + m3/4)*L*s4*th4d^2 ...
             + 2*m2*c1*yd*th1d;
    
    Phi(2) = F2 + m2*y*th1d^2 + m2*g*c1;
    
    Phi(3) = -(5*m1/128 + m3/16)*L^2*s12*th2d^2 ...
             -(3*m1/128 + m3/16)*L^2*s13*th3d^2 ...
             -(m1/128 + m3/16)*L^2*s14*th4d^2 ...
             -2*m2*y*yd*th1d - k*th1 + m2*g*s1*y ...
             + (7*m1/32 + m3/4)*g*L*s1;
    
    Phi(4) = (5*m1/128 + m3/16)*L^2*s12*th1d^2 ...
             -(3*m1/128 + m3/16)*L^2*s23*th3d^2 ...
             -(m1/128 + m3/16)*L^2*s24*th4d^2 ...
             - k*th2 + (5*m1/32 + m3/4)*g*L*s2;
    
    Phi(5) = (3*m1/128 + m3/16)*L^2*(s13*th1d^2 + s23*th2d^2) ...
             -(m1/128 + m3/16)*L^2*s34*th4d^2 ...
             - k*th3 + (3*m1/32 + m3/4)*g*L*s3;
    
    Phi(6) = (m1/128 + m3/16)*L^2*(s14*th1d^2 + s24*th2d^2 + s34*th3d^2) ...
             - k*th4 + (m1/32 + m3/4)*g*L*s4;
    
    % Solve for accelerations
    qdd = Mmat \ Phi;
    
    % Build state derivative
    dstate = zeros(12,1);
    dstate(1) = state(7);   % xd
    dstate(2) = state(8);   % yd
    dstate(3) = state(9);   % th1d
    dstate(4) = state(10);  % th2d
    dstate(5) = state(11);  % th3d
    dstate(6) = state(12);  % th4d
    dstate(7) = qdd(1);     % xdd
    dstate(8) = qdd(2);     % ydd
    dstate(9) = qdd(3);     % th1dd
    dstate(10) = qdd(4);    % th2dd
    dstate(11) = qdd(5);    % th3dd
    dstate(12) = qdd(6);    % th4dd
end

%% Second Pendulum Dynamics (L/4 <= y < L/2)
function dstate = second_pendulum(~, state, P, Forces)
    % Extract parameters
    L = P.L; M = P.M; m1 = P.m1; m2 = P.m2; m3 = P.m3; g = P.g; k = P.k;
    
    % Extract state variables
    y = state(2);
    th1 = state(3); th2 = state(4); th3 = state(5); th4 = state(6);
    yd = state(8);
    th1d = state(9); th2d = state(10); th3d = state(11); th4d = state(12);
    
    % Compute all trigonometric terms at once
    t = compute_trig(th1, th2, th3, th4);
    s1 = t.s1; c1 = t.c1; s2 = t.s2; c2 = t.c2;
    s3 = t.s3; c3 = t.c3; s4 = t.s4; c4 = t.c4;
    s12 = t.s12; c12 = t.c12; s13 = t.s13; c13 = t.c13;
    s14 = t.s14; c14 = t.c14; s23 = t.s23; c23 = t.c23;
    s24 = t.s24; c24 = t.c24; s34 = t.s34; c34 = t.c34;
    
    % Forces
    F1 = Forces(1);
    F2 = Forces(2);
    
    % Coefficients
    k1 = M + m1 + m2 + m3;
    k2 = (7*m1/32 + (m2 + m3)/4)*L;
    k3 = 5*L*m1/32 - (m2 - m3)*L/4 + m2*y;
    k4 = (3*m1/32 + m3/4)*L;
    k5 = (m1/32 + m3/4)*L;
    k6 = (5*m1*L^2/128 + (m2*L^2 + m3*L^2)/16 + m2*L*y/4);
    k7 = (13*m1/128 + m3/16)*L^2;
    k8 = 5*m1*L^2/128 - (m2 - m3)*L^2/16 + m2*L*y/4;
    k9 = (3*m1/128 + m3/16)*L^2;
    k10 = (m1/128 + m3/16)*L^2;
    
    % Build mass matrix
    Mmat = zeros(6,6);
    Mmat(1,1) = k1;
    Mmat(1,2) = m2*s2;
    Mmat(1,3) = k2*c1;
    Mmat(1,4) = k3*c2;
    Mmat(1,5) = k4*c3;
    Mmat(1,6) = k5*c4;
    
    Mmat(2,1) = m2*s2;
    Mmat(2,2) = m2;
    Mmat(2,3) = -m2*L/4*s12;
    
    Mmat(3,1) = k2*c1;
    Mmat(3,2) = -m2*L/4*s12;
    Mmat(3,3) = (13*m1/256+(m2+m3)/16)*L^2;
    Mmat(3,4) = k6*c12;
    Mmat(3,5) = k7*c13;
    Mmat(3,6) = k7*c14;
    
    Mmat(4,1) = k3*c2;
    Mmat(4,3) = k8*c12;
    Mmat(4,4) = 9*m1*L^2/256 + m2*L^2/16 + m3*L^2/16 + m2*y^2 - m2*L*y/2;
    Mmat(4,5) = k9*c23;
    Mmat(4,6) = k10*c24;
    
    Mmat(5,1) = k4*c3;
    Mmat(5,3) = k9*c13;
    Mmat(5,4) = k9*c23;
    Mmat(5,5) = (5*m1/256 + m3/16)*L^2;
    Mmat(5,6) = k10*c34;
    
    Mmat(6,1) = k5*c4;
    Mmat(6,3) = k10*c14;
    Mmat(6,4) = k10*c24;
    Mmat(6,5) = k10*c34;
    Mmat(6,6) = (m1/256 + m3/16)*L^2;
    
    % Build forcing vector
    Phi = zeros(6,1);
    Phi(1) = F1 + k2*s1*th1d^2 + k3*s2*th2d^2 - 2*m2*c2*yd*th2d ...
             + k4*s3*th3d^2 + k5*s4*th4d^2;
    
    Phi(2) = F2 + m2*L/4*c12*th1d^2 + m2*g*c2 + (-m2*L/4 + m2*y)*th2d^2;
    
    Phi(3) = -k6*s12*th2d^2 - m2*L/2*c12*yd*th2d - k7*(s13*th3d^2 + s14*th4d^2) ...
             - k*th1 + (7*m1/32+(m2+m3)/4)*g*L*s1;
    
    Phi(4) = k8*s12*th1d^2 - (2*m2*y - m2*L/2)*yd*th2d - k9*s23*th3d^2 ...
             - k10*s24*th4d^2 - k*th2 + (m2*y+5*m1*L/32-(m2-m3)*L/4)*g*s2;
    
    Phi(5) = (3*m1/128+m3/16)*L^2*(s13*th1d^2+s23*th2d^2) ...
             - (m1/128+m3/16)*L^2*s34*th4d^2 - k*th3 + (3*m1/32+m3/4)*g*L*s3;
    
    Phi(6) = (m1/128+m3/16)*L^2*(s14*th1d^2+s24*th2d^2+s34*th3d^2) ...
             - k*th4 + (m1/32+m3/4)*g*L*s4;
    
    % Solve for accelerations
    qdd = Mmat \ Phi;
    
    % Build state derivative
    dstate = zeros(12,1);
    dstate(1) = state(7);
    dstate(2) = state(8);
    dstate(3) = state(9);
    dstate(4) = state(10);
    dstate(5) = state(11);
    dstate(6) = state(12);
    dstate(7) = qdd(1);
    dstate(8) = qdd(2);
    dstate(9) = qdd(3);
    dstate(10) = qdd(4);
    dstate(11) = qdd(5);
    dstate(12) = qdd(6);
end

%% Third Pendulum Dynamics (L/2 <= y < 3L/4)
function dstate = third_pendulum(~, state, P, Forces)
    % Extract parameters
    L = P.L; M = P.M; m1 = P.m1; m2 = P.m2; m3 = P.m3; g = P.g; k = P.k;
    
    % Extract state variables
    y = state(2);
    th1 = state(3); th2 = state(4); th3 = state(5); th4 = state(6);
    yd = state(8);
    th1d = state(9); th2d = state(10); th3d = state(11); th4d = state(12);
    
    % Compute all trigonometric terms at once
    t = compute_trig(th1, th2, th3, th4);
    s1 = t.s1; c1 = t.c1; s2 = t.s2; c2 = t.c2;
    s3 = t.s3; c3 = t.c3; s4 = t.s4; c4 = t.c4;
    s12 = t.s12; c12 = t.c12; s13 = t.s13; c13 = t.c13;
    s14 = t.s14; c14 = t.c14; s23 = t.s23; c23 = t.c23;
    s24 = t.s24; c24 = t.c24; s34 = t.s34; c34 = t.c34;
    
    % Forces
    F1 = Forces(1);
    F2 = Forces(2);
    
    % Coefficients
    k1 = M + m1 + m2 + m3;
    k2 = (7*m1/32 + (m2 + m3)/4)*L;
    k3 = (5*m1/32 + (m2 + m3)/4)*L;
    k4 = 3*m1*L/32 - m2*L/2 + m3*L/4 + m2*y;
    k5 = (3*m1/32 + m3/4)*L;
    k6 = m2*L/4;
    k7 = (5*m1/128 + (m2 + m3)/16)*L^2;
    k8 = 3*m1*L^2/128 - m2*L^2/8 + m3*L^2/16 + m2*L*y/4;
    k9 = (m1/128 + m3/16)*L^2;
    k10 = 5*m1/128 + m2/16 + m3/16;
    k11 = 3*m1*L^2/128 - m2*L^2/8 + m3*L^2/16 + m2*L*y/4;
    k12 = 3*m1*L^2/128 - m2*L^2/8 + m3*L^2/16 + m3*L*y/4;
    
    % Build mass matrix
    Mmat = zeros(6,6);
    Mmat(1,1) = k1;
    Mmat(1,2) = m2*s3;
    Mmat(1,3) = k2*c1;
    Mmat(1,4) = k3*c2;
    Mmat(1,5) = k4*c3;
    Mmat(1,6) = k5*c4;
    
    Mmat(2,1) = m2*s3;
    Mmat(2,2) = m2;
    Mmat(2,3) = -m2*L*s13/4;
    Mmat(2,4) = -m2*L*s23/4;
    
    Mmat(3,1) = k2*c1;
    Mmat(3,2) = -m2*L*s13/4;
    Mmat(3,3) = (13*m1/256 + m2/16 + m3/16)*L^2;
    Mmat(3,4) = k7*c12;
    Mmat(3,5) = k12*c13;
    Mmat(3,6) = k9*c14;
    
    Mmat(4,1) = k3*c2;
    Mmat(4,2) = -m2*L*s23/4;
    Mmat(4,3) = k10*c12;
    Mmat(4,4) = (9*m1/256 + (m2 + m3)/16)*L^2;
    Mmat(4,5) = k11*c23;
    Mmat(4,6) = k9*c24;
    
    Mmat(5,1) = (m2*y + 3*m1*L/32 - m2*L/2 + m3*L/4)*c3;
    Mmat(5,3) = k8*c13;
    Mmat(5,4) = k8*c23;
    Mmat(5,5) = m2*y^2 + 5*m1*L^2/256 + m2*L^2/4 + m3*L^2/16 - m2*L*y;
    Mmat(5,6) = k9*c34;
    
    Mmat(6,1) = (m1/32 + m3/4)*L*c4;
    Mmat(6,3) = (m1/128 + m3/16)*L^2*c14;
    Mmat(6,4) = (m1/128 + m3/16)*L^2*c24;
    Mmat(6,5) = (m1/128 + m3/16)*L^2*c34;
    Mmat(6,6) = (m1/256 + m3/16)*L^2;
    
    % Build forcing vector
    Phi = zeros(6,1);
    Phi(1) = F1 + k2*s1*th1d^2 + k3*s2*th2d^2 + k4*s3*th3d^2 ...
             - 2*m2*c3*yd*th3d + k5*s4*th4d^2;
    
    Phi(2) = F2 + k6*(c13*th1d^2 + c23*th2d^2) - (m2*L/2 - m2*y)*th3d^2 + m2*g*c3;
    
    Phi(3) = -k7*s12*th2d^2 - k12*s13*th3d^2 - k9*s14*th4d^2 ...
             - k*th1 - m2*L/2*c13*yd*th3d + (7*m1/32 + (m2 + m3)/4)*g*L*s1;
    
    Phi(4) = k10*s12*th1d^2 - k*th2 - k11*s23*th3d^2 ...
             - (m1/128 + m3/16)*L^2*s24*th4d^2 - m2*L/2*c23*yd*th3d + k3*g*s2;
    
    Phi(5) = k8*(s13*th1d^2 + s23*th2d^2) - k9*s34*th4d^2 - k*th3 ...
             - (2*m2*y - m2*L)*yd*th3d ...
             + (m2*y + 3*m1*L/32 - m2*L/2 + m3*L/4)*g*s3;
    
    Phi(6) = (m1/128 + m3/16)*L^2*(s14*th1d^2 + s24*th2d^2 + s34*th3d^2) ...
             - k*th4 + (m1/32 + m3/4)*g*L*s4;
    
    % Solve for accelerations
    qdd = Mmat \ Phi;
    
    % Build state derivative
    dstate = zeros(12,1);
    dstate(1) = state(7);
    dstate(2) = state(8);
    dstate(3) = state(9);
    dstate(4) = state(10);
    dstate(5) = state(11);
    dstate(6) = state(12);
    dstate(7) = qdd(1);
    dstate(8) = qdd(2);
    dstate(9) = qdd(3);
    dstate(10) = qdd(4);
    dstate(11) = qdd(5);
    dstate(12) = qdd(6);
end

%% Fourth Pendulum Dynamics (3L/4 <= y <= L)
function dstate = fourth_pendulum(~, state, P, Forces)
    % Extract parameters
    L = P.L; M = P.M; m1 = P.m1; m2 = P.m2; m3 = P.m3; g = P.g; k = P.k;
    
    % Extract state variables
    y = state(2);
    th1 = state(3); th2 = state(4); th3 = state(5); th4 = state(6);
    yd = state(8);
    th1d = state(9); th2d = state(10); th3d = state(11); th4d = state(12);
    
    % Compute all trigonometric terms at once
    t = compute_trig(th1, th2, th3, th4);
    s1 = t.s1; c1 = t.c1; s2 = t.s2; c2 = t.c2;
    s3 = t.s3; c3 = t.c3; s4 = t.s4; c4 = t.c4;
    s12 = t.s12; c12 = t.c12; s13 = t.s13; c13 = t.c13;
    s14 = t.s14; c14 = t.c14; s23 = t.s23; c23 = t.c23;
    s24 = t.s24; c24 = t.c24; s34 = t.s34; c34 = t.c34;
    
    % Forces
    F1 = Forces(1);
    F2 = Forces(2);
    
    % Coefficients
    k1 = M + m1 + m2 + m3;
    k2 = (7*m1/32 + (m2 + m3)/4)*L;
    k3 = (5*m1/32 + (m2 + m3)/4)*L;
    k4 = (3*m1/32 + (m2 + m3)/4)*L;
    k5 = m1*L/32 - 3*m2*L/4 + m3*L/4 + m2*y;
    k6 = (5*m1/128 + (m2 + m3)/16)*L^2;
    k7 = (3*m1/128 + (m2 + m3)/16)*L^2;
    k8 = 13*m1*L^2/128 - 3*m2*L^2/16 + m3*L^2/16 + m2*y*L/4;
    k9 = m1*L^2/128 - 3*m2*L^2/16 + m3*L^2/16 + m2*y*L/4;
    k10 = m1*L^2/128 - (3*m2 + m3)*L^2/16 + m2*y*L/4;
    
    % Build mass matrix
    Mmat = zeros(6,6);
    Mmat(1,1) = k1;
    Mmat(1,2) = m2*s4;
    Mmat(1,3) = k2*c1;
    Mmat(1,4) = k3*c2;
    Mmat(1,5) = k4*c3;
    Mmat(1,6) = k5*c4;
    
    Mmat(2,1) = m2*s4;
    Mmat(2,2) = m2;
    Mmat(2,3) = -m2*L*s14/4;
    Mmat(2,4) = -m2*L*s24/4;
    Mmat(2,5) = -m2*L*s34/4;
    
    Mmat(3,1) = k2*c1;
    Mmat(3,2) = -m2*L*s14/4;
    Mmat(3,3) = (13*m1/256 + m2/16 + m3/16)*L^2;
    Mmat(3,4) = k6*c12;
    Mmat(3,5) = k7*c13;
    Mmat(3,6) = k8*c14;
    
    Mmat(4,1) = k3*c2;
    Mmat(4,2) = -m2*L*s24/4;
    Mmat(4,3) = k6*c12;
    Mmat(4,4) = (9*m1/256 + (m2 + m3)/16)*L^2;
    Mmat(4,5) = k7*c23;
    Mmat(4,6) = k9*c24;
    
    Mmat(5,1) = k4*c3;
    Mmat(5,2) = -m2*L*s34/4;
    Mmat(5,3) = k7*c13;
    Mmat(5,4) = k7*c23;
    Mmat(5,5) = (5*m1/256 + (m2 + m3)/16)*L^2;
    Mmat(5,6) = k10*c34;
    
    Mmat(6,1) = (m2*y + m1*L/32 - 3*m2*L/4 + m3*L/4)*c4;
    Mmat(6,3) = k9*c14;
    Mmat(6,4) = k9*c24;
    Mmat(6,5) = k9*c34;
    Mmat(6,6) = m2*y^2 + m1*L^2/256 + (9*m2 + m3)*L^2/16 - 3*m2*y*L/2;
    
    % Build forcing vector
    Phi = zeros(6,1);
    Phi(1) = F1 + k2*s1*th1d^2 + k3*s2*th2d^2 + k4*s3*th3d^2 ...
             + k5*s4*th4d^2 - 2*m2*c4*yd*th4d;
    
    Phi(2) = F2 + m2*L/4*(c14*th1d^2 + c24*th2d^2 + c34*th3d^2) ...
             - (3*m2*L/4 - m2*y)*th4d^2 + m2*g*c4;
    
    Phi(3) = -k6*s12*th2d^2 - k7*s13*th3d^2 - k8*s14*th4d^2 ...
             - m2*L/2*c14*yd*th4d - k*th1 + (7*m1/32 + (m2 + m3)/4)*g*L*s1;
    
    Phi(4) = (5*m1/128 + (m2 + m3)/16)*s12*th1d^2 - k7*s23*th3d^2 ...
             - m2*L*c24*yd*th4d + k*th1 - k9*s24*th4d^2 - k*th2 + k3*g*s2;
    
    Phi(5) = k7*(s13*th1d^2 + s23*th2d^2) - k10*s34*th4d^2 ...
             - m2*L/2*c34*yd*th4d - k*th3 + k4*g*s3;
    
    Phi(6) = -(2*m2*y - 3*m2*L/2)*yd*th4d ...
             + k9*(s14*th1d^2 + s24*th2d^2 + s34*th3d^2) ...
             - k*th4 + (m2*y + m1*L/32 - 3*m2*L/4 + m3*L/4)*g*s4;
    
    % Solve for accelerations
    qdd = Mmat \ Phi;
    
    % Build state derivative
    dstate = zeros(12,1);
    dstate(1) = state(7);
    dstate(2) = state(8);
    dstate(3) = state(9);
    dstate(4) = state(10);
    dstate(5) = state(11);
    dstate(6) = state(12);
    dstate(7) = qdd(1);
    dstate(8) = qdd(2);
    dstate(9) = qdd(3);
    dstate(10) = qdd(4);
    dstate(11) = qdd(5);
    dstate(12) = qdd(6);
end