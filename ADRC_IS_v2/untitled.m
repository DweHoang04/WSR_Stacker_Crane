syms t
u = heaviside(t);
fplot(u, [-5, 5]);
xlabel('t'); ylabel('u(t)');
title('Unit Step Function using heaviside');
grid on;
