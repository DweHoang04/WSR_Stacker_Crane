close all;
B = [100, 100, 600 300];
figure('Position', B);
figure(1)
hold on;
grid on
plot(out.tout,out.F1,'linewidth',1.2);
plot(out.tout,out.F11,'-.','color',[0.12,0.66,0.12],'linewidth',1.2);
plot(out.tout,out.F12,'color',[0.85,0.33,0.10],'linewidth',1.2);
plot(out.tout,out.F13,'--','color',[0,0,0],'linewidth',1.2);
legend('PID', 'PID+ZV','PID+ZVD','PID+ETM4');
xlabel('Thời gian (s)');
ylabel('F1 (N)');
axis([0 10 -5 350]);

% do thi vi tri
figure('Position', B);
figure(2)
hold on 
grid on 
plot(out.tout,out.X1,'linewidth',1.2);
plot(out.tout,out.X11,'-.','color',[0.12,0.66,0.12],'linewidth',1.2);
plot(out.tout,out.X12,'color',[0.85,0.33,0.10],'linewidth',1.2);
plot(out.tout,out.X13,'--k','linewidth',1.5);
xlabel('Thời gian (s)');
ylabel('x1 (m)');
legend('PID', 'PID+ZV','PID+ZVD','PID+ETM4');
axis([0 10 0 0.25]);


% do thi dao dong diem giua
figure('Position', [100 100 700 400]);
hold on 
grid on
subplot(2,2,1)
plot(out.tout,out.W_MID,'linewidth',1.2);
legend('PID');
xlabel('Thời gian (s)');
ylabel('ω(L/2,t)(m)')
axis([0 10 -0.03 0.03]);
grid on

subplot(2,2,2)
plot(out.tout,out.W_MID1,'color',[0.12,0.66,0.12],'linewidth',1);
legend('PID+ZV');
xlabel('Thời gian (s)');
ylabel('ω(L/2,t)(m)')
axis([0 10 -0.03 0.03]);
grid on

subplot(2,2,3)
plot(out.tout,out.W_MID2,'color',[0.85,0.33,0.10],'linewidth',1);
legend('PID+ZVD');
xlabel('Thời gian (s)');
ylabel('ω(L/2,t)(m)')
axis([0 10 -0.03 0.03]);
grid on

subplot(2,2,4)
plot(out.tout,out.W_MID3,'k','linewidth',1.2);
legend('PID+ETM4');
xlabel('Thời gian (s)');
ylabel('ω(L/2,t)(m)')
axis([0 10 -0.03 0.03]);
grid on

% plot(out.tout,out.W_MID1,'--r','linewidth',1.2);
% plot(out.tout,out.W_MID2,'-.k','linewidth',1.2);
% plot(out.tout,out.W_MID3,'-m','linewidth',1.2);
% xlabel('Thời gian (s)');
% ylabel('ω(L/2, t)(m)');
% legend('PID', 'PID+ZV','PID+ZVD','PID+ETM4');

% do thi dao dong diem cuoi
figure('Position', [100 100 700 400]);
hold on 
grid on
subplot(2,2,1)
plot(out.tout,out.WL,'linewidth',1.2);
legend('PID');
xlabel('Thời gian (s)');
ylabel('ω(L,t)(m)')
axis([0 10 -0.05 0.05]);
grid on

subplot(2,2,2)
plot(out.tout,out.WL1,'color',[0.12,0.66,0.12],'linewidth',1);
legend('PID+ZV');
xlabel('Thời gian (s)');
ylabel('ω(L,t)(m)')
axis([0 10 -0.05 0.05]);
grid on

subplot(2,2,3)
plot(out.tout,out.WL2,'color',[0.85,0.33,0.10],'linewidth',1);
legend('PID+ZVD');
xlabel('Thời gian (s)');
ylabel('ω(L,t)(m)')
axis([0 10 -0.05 0.05]);
grid on

subplot(2,2,4)
plot(out.tout,out.WL3,'k','linewidth',1.2);
legend('PID+ETM4');
xlabel('Thời gian (s)');
ylabel('ω(L,t)(m)')
axis([0 10 -0.05 0.05]);
grid on


% figure('Position', B);
% %figure(4)
% hold on 
% grid on 
% plot(out.tout,out.WL,'linewidth',0.7);
% plot(out.tout,out.WL1,'linewidth',0.7);
% plot(out.tout,out.WL2,'linewidth',0.7);
% plot(out.tout,out.WL3,'linewidth',0.7);
% xlabel('Thời gian (s)');
% ylabel('ω(L,t) (m)');
% legend('Middle point','End point');
% axis([0 15 -0.04 0.04])
% legend('PID', 'PID+ZV','PID+ZVD','PID+ETM4');

%do thi pho tan so
% figure('Position', B);
% tmax = 15;
% r = 5000;
% dt = tmax/(r);
% Fs=1/dt;
% time=0:dt:(r-2)*dt-dt;
% l1=r;
% grid on
% hold on
% 
% fft_w=fft(out.W_MID,l1)*(2/l1);
% abs_w=abs(fft_w);
% freq=0:(1/time(end)):Fs/2-(1/time(end));
% plot(freq,abs_w(1:length(freq)),LineWidth=1.2); 
% 
% fft_w=fft(out.WL,l1)*(2/l1);
% abs_w=abs(fft_w);
% freq=0:(1/time(end)):Fs/2-(1/time(end));
% plot(freq,abs_w(1:length(freq)),LineWidth=1.2); 
% xlabel('Tần số (Hz)');
% ylabel('Biên độ (m)')
% legend('Middle point','End point');
% A = [0 10 0 0.02 ];
% axis(A);

