%% Stacker Crane Simulation
% Based on: "Robust Command Shaped Vibration Control for Stacker Crane"
% by Li et al., IEEE Transactions on Industrial Electronics, 2024

clear; clc; close all;

%% Initial Conditions
Parameters = struct();
Parameters.L = 0.4;     % Total mast length (m)
Parameters.M = 1.5;     % Trolley mass (kg)
Parameters.m1 = 0.075;  % First pendulum mass (kg)
Parameters.m2 = 0.15;   % Payload mass (kg)
Parameters.m3 = 0.03;   % Tip mass (kg)
Parameters.g = 9.81;    % Gravitational acceleration (m/s^2)
Parameters.k = 3.26;    % Virtual spring stiffness (N*m/rad)

%% Initial Conditions
% State vector: [x, y, th1, th2, th3, th4, xd, yd, th1d, th2d, th3d, th4d]
Y0 = [0;        % x:    Initial trolley position (m)
          0.05;     % y:    Initial payload height (m)
          0;        % th1:  Initial angle of pendulum 1 (rad)
          0;        % th2:  Initial angle of pendulum 2 (rad)
          0;        % th3:  Initial angle of pendulum 3 (rad)
          0;        % th4:  Initial angle of pendulum 4 (rad)
          0;        % xd:   Initial trolley velocity (m/s)
          0;        % yd:   Initial payload velocity (m/s)
          0;        % th1d: Initial angular velocity 1 (rad/s)
          0;        % th2d: Initial angular velocity 2 (rad/s)
          0;        % th3d: Initial angular velocity 3 (rad/s)
          0];       % th4d: Initial angular velocity 4 (rad/s)

%% Simulation Parameters
tspan = [0 10];         % Simulation time span
tstart = tspan(1);      % Starting time
tend = tspan(end);      % Ending time

%% Setup starting equation
t = tstart;
Y = Y0;
fcn = @first_pendulum;
opt = odeset('Events', @InFirstPend);

%% Main while loop
while (t(end) < tend)
    % Run integration until event functions stop it:
    [at, aY, ate, aye, aie] = ode45(fcn, [t(end), tend], Y(end,:), opt);

    % Append the new trajectory:
    t = cat(1, t, at(2:end));
    Y = cat(1, aY(2:end, :));

    if 
end

%% ========================================================================
%  SEPARATE DYNAMIC EQUATIONS FOR EACH PENDULUM POSITION
%  ========================================================================

%% First Pendulum Dynamics (0 <= y < L/4)
function dY = first_pendulum(t, Y, P, Forces)
    % Extract parameters
    L = P.L; M = P.M; m1 = P.m1; m2 = P.m2; m3 = P.m3; g = P.g; k = P.k;
    
    % Extract Y variables
    y = Y(2);
    th1 = Y(3); th2 = Y(4); th3 = Y(5); th4 = Y(6);
    yd = Y(8);
    th1d = Y(9); th2d = Y(10); th3d = Y(11); th4d = Y(12);
    
    % Trigonometric functions
    s1 = sin(th1); c1 = cos(th1);
    s2 = sin(th2); c2 = cos(th2);
    s3 = sin(th3); c3 = cos(th3);
    s4 = sin(th4); c4 = cos(th4);
    
    s12 = sin(th1-th2); c12 = cos(th1-th2);
    s13 = sin(th1-th3); c13 = cos(th1-th3);
    s14 = sin(th1-th4); c14 = cos(th1-th4);
    s23 = sin(th2-th3); c23 = cos(th2-th3);
    s24 = sin(th2-th4); c24 = cos(th2-th4);
    s34 = sin(th3-th4); c34 = cos(th3-th4);
    
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
    
    % Build Y derivative
    dY = zeros(12,1);
    dY(1) = Y(7);   % xd
    dY(2) = Y(8);   % yd
    dY(3) = Y(9);   % th1d
    dY(4) = Y(10);  % th2d
    dY(5) = Y(11);  % th3d
    dY(6) = Y(12);  % th4d
    dY(7) = qdd(1);     % xdd
    dY(8) = qdd(2);     % ydd
    dY(9) = qdd(3);     % th1dd
    dY(10) = qdd(4);    % th2dd
    dY(11) = qdd(5);    % th3dd
    dY(12) = qdd(6);    % th4dd
end

%% Second Pendulum Dynamics (L/4 <= y < L/2)
function dY = second_pendulum(t, Y, P, Forces)
    % Extract parameters
    L = P.L; M = P.M; m1 = P.m1; m2 = P.m2; m3 = P.m3; g = P.g; k = P.k;
    
    % Extract Y variables
    y = Y(2);
    th1 = Y(3); th2 = Y(4); th3 = Y(5); th4 = Y(6);
    yd = Y(8);
    th1d = Y(9); th2d = Y(10); th3d = Y(11); th4d = Y(12);
    
    % Trigonometric functions
    s1 = sin(th1); c1 = cos(th1);
    s2 = sin(th2); c2 = cos(th2);
    s3 = sin(th3); c3 = cos(th3);
    s4 = sin(th4); c4 = cos(th4);
    
    s12 = sin(th1-th2); c12 = cos(th1-th2);
    s13 = sin(th1-th3); c13 = cos(th1-th3);
    s14 = sin(th1-th4); c14 = cos(th1-th4);
    s23 = sin(th2-th3); c23 = cos(th2-th3);
    s24 = sin(th2-th4); c24 = cos(th2-th4);
    s34 = sin(th3-th4); c34 = cos(th3-th4);
    
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
    
    % Build Y derivative
    dY = zeros(12,1);
    dY(1) = Y(7);   % xd
    dY(2) = Y(8);   % yd
    dY(3) = Y(9);   % th1d
    dY(4) = Y(10);  % th2d
    dY(5) = Y(11);  % th3d
    dY(6) = Y(12);  % th4d
    dY(7) = qdd(1);     % xdd
    dY(8) = qdd(2);     % ydd
    dY(9) = qdd(3);     % th1dd
    dY(10) = qdd(4);    % th2dd
    dY(11) = qdd(5);    % th3dd
    dY(12) = qdd(6);    % th4dd
end

%% Third Pendulum Dynamics (L/2 <= y < 3L/4)
function dY = third_pendulum(t, Y, P, Forces)
    % Extract parameters
    L = P.L; M = P.M; m1 = P.m1; m2 = P.m2; m3 = P.m3; g = P.g; k = P.k;
    
    % Extract Y variables
    y = Y(2);
    th1 = Y(3); th2 = Y(4); th3 = Y(5); th4 = Y(6);
    yd = Y(8);
    th1d = Y(9); th2d = Y(10); th3d = Y(11); th4d = Y(12);
    
    % Trigonometric functions
    s1 = sin(th1); c1 = cos(th1);
    s2 = sin(th2); c2 = cos(th2);
    s3 = sin(th3); c3 = cos(th3);
    s4 = sin(th4); c4 = cos(th4);
    
    s12 = sin(th1-th2); c12 = cos(th1-th2);
    s13 = sin(th1-th3); c13 = cos(th1-th3);
    s14 = sin(th1-th4); c14 = cos(th1-th4);
    s23 = sin(th2-th3); c23 = cos(th2-th3);
    s24 = sin(th2-th4); c24 = cos(th2-th4);
    s34 = sin(th3-th4); c34 = cos(th3-th4);
    
    % Forces
    F1 = Forces(1);
    F2 = Forces(2);
    
    % Coefficients
    k1  = M + m1 + m2 + m3;
    k2  = (7*m1/32 + (m2 + m3)/4)*L;
    k3  = (5*m1/32 + (m2 + m3)/4)*L;
    k4  = 3*m1*L/32 - m2*L/2 + m3*L/4 + m2*y;
    k5  = (m1/32 + m3/4)*L; % Corrected from original code (was k5 for 3rd pendulum)
    k6  = m2*L/4;
    k7  = (5*m1/128 + (m2 + m3)/16)*L^2;
    k8  = 3*m1*L^2/128 - m2*L^2/8 + m3*L^2/16 + m2*L*y/4;
    k9  = (m1/128 + m3/16)*L^2;
    k10 = (5*m1/128 + (m2+m3)/16)*L^2; % Corrected from original code
    k11 = 3*m1*L^2/128 - m2*L^2/8 + m3*L^2/16 + m2*L*y/4;
    k12 = (3*m1/128 + (m2+m3)/16)*L^2; % Corrected from original code
    
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
    Mmat(3,3) = (13*m1/256 + (m2 + m3)/16)*L^2;
    Mmat(3,4) = k7*c12;
    Mmat(3,5) = k12*c13;
    Mmat(3,6) = k9*c14;
    
    Mmat(4,1) = k3*c2;
    Mmat(4,2) = -m2*L*s23/4;
    Mmat(4,3) = k7*c12; % k10 in original code was k7
    Mmat(4,4) = (9*m1/256 + (m2 + m3)/16)*L^2;
    Mmat(4,5) = k11*c23;
    Mmat(4,6) = k9*c24;
    
    Mmat(5,1) = (m2*y + 3*m1*L/32 - m2*L/2 + m3*L/4)*c3;
    Mmat(5,3) = k8*c13;
    Mmat(5,4) = k8*c23;
    Mmat(5,5) = (m2*y^2 + 5*m1*L^2/256 + m2*L^2/4 + m3*L^2/16 - m2*L*y);
    Mmat(5,6) = k9*c34;
    
    Mmat(6,1) = (m1/32 + m3/4)*L*c4;
    Mmat(6,3) = (m1/128 + m3/16)*L^2*c14;
    Mmat(6,4) = (m1/128 + m3/16)*L^2*c24;
    Mmat(6,5) = (m1/128 + m3/16)*L^2*c34;
    Mmat(6,6) = (m1/256 + m3/16)*L^2;

    % Build forcing vector
    Phi = zeros(6,1);
    Phi(1) = F1 + k2*s1*th1d^2 + k3*s2*th2d^2 + k4*s3*th3d^2 - 2*m2*c3*yd*th3d + k5*s4*th4d^2;
    Phi(2) = F2 + k6*(c13*th1d^2 + c23*th2d^2) - (m2*L/2 - m2*y)*th3d^2 + m2*g*c3;
    Phi(3) = -k7*s12*th2d^2 - k12*s13*th3d^2 - k9*s14*th4d^2 - k*th1 - m2*L/2*c13*yd*th3d + (7*m1/32 + (m2 + m3)/4)*g*L*s1;
    Phi(4) = k7*s12*th1d^2 - k*th2 - k11*s23*th3d^2 - (m1/128 + m3/16)*L^2*s24*th4d^2 - m2*L/2*c23*yd*th3d + k3*g*s2;
    Phi(5) = k8*(s13*th1d^2 + s23*th2d^2) - k9*s34*th4d^2 - k*th3 - (2*m2*y - m2*L)*yd*th3d + (m2*y + 3*m1*L/32 - m2*L/2 + m3*L/4)*g*s3;
    Phi(6) = (m1/128 + m3/16)*L^2*(s14*th1d^2 + s24*th2d^2 + s34*th3d^2) - k*th4 + (m1/32 + m3/4)*g*L*s4;

    % Solve for accelerations
    qdd = Mmat \ Phi;
    
    % Build Y derivative
    dY = zeros(12,1);
    dY(1) = Y(7);   % xd
    dY(2) = Y(8);   % yd
    dY(3) = Y(9);   % th1d
    dY(4) = Y(10);  % th2d
    dY(5) = Y(11);  % th3d
    dY(6) = Y(12);  % th4d
    dY(7) = qdd(1);     % xdd
    dY(8) = qdd(2);     % ydd
    dY(9) = qdd(3);     % th1dd
    dY(10) = qdd(4);    % th2dd
    dY(11) = qdd(5);    % th3dd
    dY(12) = qdd(6);    % th4dd
end

%% Fourth Pendulum Dynamics (3L/4 <= y <= L)
function dY = fourth_pendulum(t, Y, P, Forces)
    % Extract parameters
    L = P.L; M = P.M; m1 = P.m1; m2 = P.m2; m3 = P.m3; g = P.g; k = P.k;
    
    % Extract Y variables
    y = Y(2);
    th1 = Y(3); th2 = Y(4); th3 = Y(5); th4 = Y(6);
    yd = Y(8);
    th1d = Y(9); th2d = Y(10); th3d = Y(11); th4d = Y(12);
    
    % Trigonometric functions
    s1 = sin(th1); c1 = cos(th1);
    s2 = sin(th2); c2 = cos(th2);
    s3 = sin(th3); c3 = cos(th3);
    s4 = sin(th4); c4 = cos(th4);
    
    s12 = sin(th1-th2); c12 = cos(th1-th2);
    s13 = sin(th1-th3); c13 = cos(th1-th3);
    s14 = sin(th1-th4); c14 = cos(th1-th4);
    s23 = sin(th2-th3); c23 = cos(th2-th3);
    s24 = sin(th2-th4); c24 = cos(th2-th4);
    s34 = sin(th3-th4); c34 = cos(th3-th4);
    
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
    Mmat(6,6) = (m2*y^2 + m1*L^2/256 + (9*m2 + m3)*L^2/16 - 3*m2*L*y/2);
    
    % Build forcing vector
    Phi = zeros(6,1);
    Phi(1) = F1 + k2*s1*th1d^2 + k3*s2*th2d^2 + k4*s3*th3d^2 + k5*s4*th4d^2 - 2*m2*c4*yd*th4d;
    Phi(2) = F2 + m2*L/4*(c14*th1d^2 + c24*th2d^2 + c34*th3d^2) - (3*m2*L/4 - m2*y)*th4d^2 + m2*g*c4;
    Phi(3) = -k6*s12*th2d^2 - k7*s13*th3d^2 - k8*s14*th4d^2 - m2*L/2*c14*yd*th4d - k*th1 + (7*m1/32 + (m2 + m3)/4)*g*L*s1;
    Phi(4) = (5*m1/128 + (m2 + m3)/16)*L^2*s12*th1d^2 - k7*s23*th3d^2 - m2*L*c24*yd*th4d + k*th1 - k9*s24*th4d^2 - k*th2 + k3*g*s2;
    Phi(5) = k7*(s13*th1d^2 + s23*th2d^2) - k10*s34*th4d^2 - m2*L/2*c34*yd*th4d - k*th3 + k4*g*s3;
    Phi(6) = -(2*m2*y - 3*m2*L/2)*yd*th4d + k9*(s14*th1d^2 + s24*th2d^2 + s34*th3d^2) - k*th4 + (m2*y + m1*L/32 - 3*m2*L/4 + m3*L/4)*g*s4;
    
    % Solve for accelerations
    qdd = Mmat \ Phi;
    
    % Build Y derivative
    dY = zeros(12,1);
    dY(1) = Y(7);   % xd
    dY(2) = Y(8);   % yd
    dY(3) = Y(9);   % th1d
    dY(4) = Y(10);  % th2d
    dY(5) = Y(11);  % th3d
    dY(6) = Y(12);  % th4d
    dY(7) = qdd(1);     % xdd
    dY(8) = qdd(2);     % ydd
    dY(9) = qdd(3);     % th1dd
    dY(10) = qdd(4);    % th2dd
    dY(11) = qdd(5);    % th3dd
    dY(12) = qdd(6);    % th4dd
end

%% Event functions
% The payload is in the first pendulum
function [value,isterminal,direction] = InFirstPend(t,y)
    value = double((y >= 0) && (y < Parameters.L/4));
    isterminal = 1;
    direction = 1;
end

% The payload is in the second pendulum
function [value,isterminal,direction] = InSecondPend(t,y)
    value = double((y >= Parameters.L/4) && (y < Parameters.L/2));
    isterminal = 1;
    direction = 1;
end

% The payload is in the third pendulum
function [value,isterminal,direction] = InThirdPend(t,y)
    value = double((y >= Parameters.L/2) && (y < 3*Parameters.L/4));
    isterminal = 1;
    direction = 1;
end

% The payload is in the fourth pendulum
function [value,isterminal,direction] = InFourthPend(t,y)
    value = (y >= 3*Parameters.L/4) && (y <= Parameters.L);
    isterminal = 1;
    direction = 1;
end
