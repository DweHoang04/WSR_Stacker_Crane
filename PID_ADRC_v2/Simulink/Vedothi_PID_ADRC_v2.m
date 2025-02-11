close all;

% Đồ thị lực điều khiển F1
B = [100, 100, 600 300];
figure('Position', B);
hold on;
grid on
plot(out.tout,out.F11,'color', [0.07,0.62,1.00],'linewidth',1.5);
plot(out.tout,out.F1,'--r','linewidth',2);
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
legend('ADRC','PID');
axis([-0.1 6 -5 26]);
xticks([0 2 4 6]);
yticks([-4 0 24]);

% Đồ thị vị trí xe con
figure('Position', B);
hold on 
grid on 
plot(out.tout,out.X11,'color', [0.07,0.62,1.00],'linewidth',1.5);
plot(out.tout,out.X1,'--r','linewidth',2);
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
axis([6 8 0.195 0.205]);
yticks([0.195 0.2 0.2023 0.205])

% Đồ thị độ lắc điểm cuối
figure('Position', [100 100 600 300]);
hold on 
grid on
plot(out.tout,out.WL1,'color', [0.07,0.62,1.00],'linewidth',1.5);
plot(out.tout,out.WL,'--r','linewidth',2);
legend('ADRC','PID');
xlabel('Thời gian (s)');
ylabel('Dao động (m)')
axis([0 10 -0.05 0.07]);
yticks([-0.04 -0.02 0 0.02 0.04 0.06]);
xticks([0 2 4 6 8 10]);
grid on

% Đồ thị phổ tần số
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

