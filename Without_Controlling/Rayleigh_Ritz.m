% clc;
clear all;

n=7;

syms x mk a b c d L EI pA g h q mh H
f = a*(x/L)^2 + b*(x/L)^3 + c*(x/L)^4 + d*(x/L)^5 + h*(x/L)^6 + q*(x/L)^7 + g*(x/L)^8;
coef = [a b c d h q g];
M = sym(zeros(n,n));
K = sym(zeros(n,n));
for i = 1:n
    for j = 1:n
        M(i,j) = int(diff(diff(f^2,coef(i)),coef(j)),x,0,L) + mh * subs(diff(diff(f^2,coef(i)),coef(j)), x, H);
        K(i,j) = int(diff(diff(diff(diff(f,x),x)^2,coef(i)),coef(j)),x,0,L);
    end
end
K = K/2 * EI;
M = M/2 * pA + mk;
lambda = double(subs(eig(inv(M)*K), [L EI pA mk mh H], [0.63 0.754 0.29655 0.04 0.86 0.2]));
freqs = sqrt(lambda)/(2*pi)