function system_of_PDEs
close all
clear
clc

% Discretization of the simulation domain
x = linspace(0,10,100);
t = linspace(0,10,100);

% Application of the Matlab partial differential equation solver 'pdepe'
m = 0;
sol = pdepe(m,@pde_func,@pde_ics,@pde_bcs,x,t);

% The solution matrices
y1 = sol(:,:,1);
y2 = sol(:,:,2);

figure(1)
surf(x,t,y1)
title('y_1(x,t)')
xlabel('Distance x')
ylabel('Time t')

figure(2)
surf(x,t,y2)
title('y_2(x,t)')
xlabel('Distance x')
ylabel('Time t')

% Equations to solve
function [c,f,s] = pde_func(x,t,y,dydx) 
% Equations arranged for the 'pdepe' solver
c = [1; 1];
f = [0.375; 0.299].*dydx;
A = sin(x);
s = [A; -A];
end

% Initial Conditions
function u0 = pde_ics(x) 
u0 = [1; 0];
end

% Boundary Conditions
function [pl,ql,pr,qr] = pde_bcs(xl,ul,xr,ur,t) 
pl = [0; ul(2)];
ql = [1; 0];
pr = [ur(1)-1; 0];
qr = [0; 1];
end

end