clc;
close all;
clear;
%--------------------------------------------------------------------------
%             3D Euler - Bernoulli beam model without control
%--------------------------------------------------------------------------

% Parameter of the Euler - Bernoulli Beam
l = 1.5; m = 1; EI = 8; EA = 14; T = 0; rho = 0.1;

% Space and Time parameters
n = 9; r = 15000;
tmax = 15;
delta_s = l/(n - 1); % Spatial step size
delta_t = tmax/(r - 1); % Temporal step size

% Create matrixes to save state data
u = zeros(n,r);
v = zeros(n,r);
w = zeros(n,r);

% Create matrixes to save disturbance at the top of the beam
du = zeros(r,1);
dv = zeros(r,1);
dw = zeros(r,1);

% Create matrixes to save displacement data
u_3D_free = zeros(r,n);
v_3D_free = zeros(r,n);
w_3D_free = zeros(r,n);

% Initial conditions
for i = 1:n
    u(i,1) = 0.5*(i - 1)*m*delta_s/l;
    v(i,1) = 0.5*(i - 1)*m*delta_s/l;
    w(i,1) = 0.5*(i - 1)*m*delta_s/l;
end

% Boundary disturbances
for j = 1:r
    du(j) = (1 + 4*sin(0.4*(j - 1)*delta_t))/5;
    dv(j) = (1 + 4*sin(0.4*(j - 1)*delta_t))/5;
    dw(j) = (1 + 4*sin(0.4*(j - 1)*delta_t))/5;
end

u(:,2) = u(:,1);
v(:,2) = v(:,1);
w(:,2) = w(:,1);

for j = 2:(r - 1)
    for i = 3:(n - 2)
        % 
        us = (u(i,j) - u(i - 1,j))/delta_s;
        vs = (v(i,j) - v(i - 1,j))/delta_s;
        ws = (w(i,j) - w(i - 1,j))/delta_s;

        % Derivative with respect to s using Finite Difference Method
        uss = (u(i + 1,j) - 2*u(i,j) + u(i - 1,j))/delta_s^2;
        vss = (v(i + 1,j) - 2*v(i,j) + v(i - 1,j))/delta_s^2;
        wss = (w(i + 1,j) - 2*w(i,j) + w(i - 1,j))/delta_s^2;
        ussss = (u(i + 2,j) - 4*u(i + 1,j) + 6*u(i,j) - 4*u(i - 1,j) + u(i - 2,j))/delta_s^4;
        vssss = (v(i + 2,j) - 4*v(i + 1,j) + 6*v(i,j) - 4*v(i - 1,j) + v(i - 2,j))/delta_s^4;

        % PDEs model of the Euler - Bernoulli beam
        S1 = T*uss - EI*ussss + EA*(wss*us + uss*ws) + 1.5*EA*us^2*uss + 0.5*EA*(uss*vs^2 + 2*us*vs*vss); % rho.diff(u,t,2)
        S2 = T*vss - EI*vssss + EA*(wss*vs + vss*ws) + 1.5*EA*vs^2*vss + 0.5*EA*(vss*us^2 + 2*vs*us*uss); % rho.diff(v,t,2)
        S3 = EA*wss + EA*us*uss + EA*vs*vss; % rho.diff(w,t,2)

        % Applying Finite - Difference Method to S1, S2 and S3
        u(i,j + 1) = 2*u(i,j) - u(i,j - 1) + delta_t^2*S1/rho;
        v(i,j + 1) = 2*v(i,j) - v(i,j - 1) + delta_t^2*S2/rho;
        w(i,j + 1) = 2*w(i,j) - w(i,j - 1) + delta_t^2*S3/rho;
    end

    % The tranverse deformation at the middle point of the beam
    u(2, j + 1) = (u(1,j + 1) + u(3,j + 1))/2;
    v(2, j + 1) = (v(1,j + 1) + v(3,j + 1))/2;
    w(2, j + 1) = (w(1,j + 1) + w(3,j + 1))/2;

    % Derivative with respect to s at the top of the beam using Finite Difference Method
    usl = (u(n,j) - u(n - 1,j))/delta_s;
    vsl = (v(n,j) - v(n - 1,j))/delta_s;
    wsl = (w(n,j) - w(n - 1,j))/delta_s;
    usssl = (-u(n,j) + 2*u(n - 1,j) - u(n - 2,j))/delta_s^3;
    vsssl = (-v(n,j) + 2*v(n - 1,j) - v(n - 2,j))/delta_s^3;

    % Boundary condition
    u(n,j + 1) = 2*u(n,j) - u(n,j - 1) - (T*usl + 0.5*EA*usl^3 + EA*usl*wsl + 0.5*EA*usl*vsl^2 - EI*usssl - du(j))*(delta_t^2)/m;
    v(n,j + 1) = 2*v(n,j) - v(n,j - 1) - (T*vsl + 0.5*EA*vsl^3 + EA*vsl*wsl + 0.5*EA*vsl*usl^2 - EI*vsssl - dv(j))*(delta_t^2)/m;
    w(n,j + 1) = 2*w(n,j) - w(n, j - 1) - (EA*wsl + 0.5*EA*usl^2 + 0.5*EA*vsl^2 - dw(j))*(delta_t^2)/m;
    
    u(n - 1,j + 1) = (u(n,j + 1) + u(n - 2,j + 1))/2;
    v(n - 1,j + 1) = (v(n,j + 1) + v(n - 2,j + 1))/2;
    w(n - 1,j + 1) = (w(n,j + 1) + w(n - 2,j + 1))/2;

    % Storing displacement data
    u_3D_free(1 + j,:) = u(:,j)';
    v_3D_free(1 + j,:) = v(:,j)';
    w_3D_free(1 + j,:) = w(:,j)';
end

u_3D_free(1,:) = u(:,1)';
v_3D_free(1,:) = v(:,1)';
w_3D_free(1,:) = w(:,1)';

% Save data for simulation
u0 = linspace(0,1,n);
t_tr = linspace(0,tmax,r);

% Plotting
figure(1);
surf(u0,t_tr,u_3D_free); view(45,30);
title({'Deformation of the beam in U direction without control in case 1'});
ylabel('t[s]','FontSize',12); xlabel('s[m]','FontSize',12); zlabel('u(s,t)[m]','FontSize',12);

figure(2);
surf(u0,t_tr,u_3D_free); view(45,30);
title({'Deformation of the beam in V direction without control in case 1'});
ylabel('t[s]','FontSize',12); xlabel('s[m]','FontSize',12); zlabel('v(s,t)[m]','FontSize',12);

figure(3);
surf(u0,t_tr,u_3D_free); view(45,30);
title({'Deformation of the beam in W direction without control in case 1'});
ylabel('t[s]','FontSize',12); xlabel('s[m]','FontSize',12); zlabel('w(s,t)[m]','FontSize',12);