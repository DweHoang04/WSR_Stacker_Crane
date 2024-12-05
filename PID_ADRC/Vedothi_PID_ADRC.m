close all;
B = [100, 100, 600 300];
figure('Position', B);
hold on;
grid on
plot(out.tout,out.F11,'color', [0.07,0.62,1.00],'linewidth',1.5);
plot(out.tout,out.F1,'--r','linewidth',2);
% plot(out.tout,out.F12,'color',[0.85,0.33,0.10],'linewidth',1.2);
% plot(out.tout,out.F13,'--','color',[0,0,0],'linewidth',1.2);
legend('ADRC','PID');
xlabel('Thời gian (s)');
ylabel('Lực điều khiển F_1 (N)');
xticks([0 2 4 6 8 10]);
yticks([0 100 200 300]);
axis([-0.1 10 -5 350]);

figure('Position', [100 100 400 150]);
hold on;
grid on
plot(out.tout,out.F11,'color', [0.07,0.62,1.00],'linewidth',1.5);
plot(out.tout,out.F1,'--r','linewidth',2);
% plot(out.tout,out.F12,'color',[0.85,0.33,0.10],'linewidth',1.2);
% plot(out.tout,out.F13,'--','color',[0,0,0],'linewidth',1.2);
legend('ADRC','PID');
% xlabel('Thời gian (s)');
% ylabel('F1 (N)');
axis([-0.1 6 -5 26]);
xticks([0 2 4 6]);
yticks([-4 0 24]);

% do thi vi tri
figure('Position', B);
hold on 
grid on 
plot(out.tout,out.X11,'color', [0.07,0.62,1.00],'linewidth',1.5);
plot(out.tout,out.X1,'--r','linewidth',2);
% plot(out.tout,out.X12,'color',[0.85,0.33,0.10],'linewidth',1.2);
% plot(out.tout,out.X13,'--k','linewidth',1.5);
xlabel('Thời gian (s)');
ylabel('Vị trí xe con (m)');
legend('ADRC','PID');
axis([0 10 0 0.25]);
xticks([0 2 4 6 8 10]);
yticks([0 0.05 0.1 0.15 0.2 0.25]);

figure('Position', [100 100 300 120]);
hold on 
grid on 
plot(out.tout,out.X11,'color', [0.07,0.62,1.00],'linewidth',1.5);
plot(out.tout,out.X1,'--r','linewidth',2);
% plot(out.tout,out.X12,'color',[0.85,0.33,0.10],'linewidth',1.2);
% plot(out.tout,out.X13,'--k','linewidth',1.5);
% xlabel('Thời gian (s)');
% ylabel('Vị trí xe con (m)');
% legend('ADRC','PID');
axis([6 8 0.195 0.205]);
yticks([0.195 0.2 0.2023 0.205])
% xticks([0 2 4 6 8 10]);
% yticks([0 0.05 0.1 0.15 0.2 0.25]);
% do thi dao dong diem giua
% figure('Position', [100 100 700 400]);
% hold on 
% grid on
% % subplot(2,2,1)
% plot(out.tout,out.W_MID1,'linewidth',1);
% plot(out.tout,out.W_MID,'linewidth',1.2);
% legend('ADRC');
% xlabel('Thời gian (s)');
% ylabel('ω(L/2,t)(m)')
% axis([0 10 -0.04 0.04]);
% grid on

% subplot(2,2,2)
% plot(out.tout,out.W_MID1,'color',[0.12,0.66,0.12],'linewidth',1);
% legend('ADRC+ZV');
% xlabel('Thời gian (s)');
% ylabel('ω(L/2,t)(m)')
% axis([0 10 -0.02 0.02]);
% grid on
% 
% subplot(2,2,3)
% plot(out.tout,out.W_MID2,'color',[0.85,0.33,0.10],'linewidth',1);
% legend('ADRC+ZVD');
% xlabel('Thời gian (s)');
% ylabel('ω(L/2,t)(m)')
% axis([0 10 -0.02 0.02]);
% grid on
% 
% subplot(2,2,4)
% plot(out.tout,out.W_MID3,'k','linewidth',1.2);
% legend('ADRC+ETM4');
% xlabel('Thời gian (s)');
% ylabel('ω(L/2,t)(m)')
% axis([0 10 -0.02 0.02]);
% grid on

% plot(out.tout,out.W_MID1,'--r','linewidth',1.2);
% plot(out.tout,out.W_MID2,'-.k','linewidth',1.2);
% plot(out.tout,out.W_MID3,'-m','linewidth',1.2);
% xlabel('Thời gian (s)');
% ylabel('ω(L/2, t)(m)');
% legend('ADRC', 'ADRC+ZV','ADRC+ZVD','ADRC+ETM4');

% do thi dao dong diem cuoi
figure('Position', [100 100 600 300]);
hold on 
grid on
% subplot(2,2,1)
plot(out.tout,out.WL1,'color', [0.07,0.62,1.00],'linewidth',1.5);
plot(out.tout,out.WL,'--r','linewidth',2);
legend('ADRC','PID');
xlabel('Thời gian (s)');
ylabel('Dao động (m)')
axis([0 10 -0.05 0.07]);
yticks([-0.04 -0.02 0 0.02 0.04 0.06]);
xticks([0 2 4 6 8 10]);
grid on

% subplot(2,2,2)
% plot(out.tout,out.WL1,'color',[0.12,0.66,0.12],'linewidth',1);
% legend('ADRC+ZV');
% xlabel('Thời gian (s)');
% ylabel('ω(L,t)(m)')
% axis([0 10 -0.038 0.038]);
% grid on
% 
% subplot(2,2,3)
% plot(out.tout,out.WL2,'color',[0.85,0.33,0.10],'linewidth',1);
% legend('ADRC+ZVD');
% xlabel('Thời gian (s)');
% ylabel('ω(L,t)(m)')
% axis([0 10 -0.038 0.038]);
% grid on
% 
% subplot(2,2,4)
% plot(out.tout,out.WL3,'k','linewidth',1.2);
% legend('ADRC+ETM4');
% xlabel('Thời gian (s)');
% ylabel('ω(L,t)(m)')
% axis([0 10 -0.038 0.038]);
% grid on


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
% legend('ADRC', 'ADRC+ZV','ADRC+ZVD','ADRC+ETM4');

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

