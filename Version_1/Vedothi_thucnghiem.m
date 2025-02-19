close all;

%vitrixe
M = xlsread('C:\Users\ADMIN\OneDrive - Hanoi University of Science and Technology\Tailieuhoctap\Stacker_crane\DATNCN\Simulation\DATA_THUCNGHIEM\Vitriv1.xlsx');
%M = xlsread('C:\Users\Admin\Downloads\abc.xlsx');
B = [100, 100, 600 300];
figure('Position', B);

M(:,8) = M(:,8)*0.01;
M(:,11) = M(:,11)*0.01;
M(:,10) = M(:,10)*0.01;
M(:,9) = M(:,9)*0.01;

hold on 
grid on 
figure(1)
% plot(out.tout,out.X1,'--r','linewidth',1.2);
plot(M(:,1),M(:,8),'color',[0.00,0.45,0.74],'linewidth',1.2);

% plot(M(:,1),M(:,9),'-.','color',[0.12,0.66,0.12],'linewidth',1.2);
% plot(M(:,1),M(:,10),'color',[0.85,0.33,0.10],'linewidth',1.2);
% plot(M(:,1),M(:,11),'--','color',[0,0,0],'linewidth',1.2);
xlabel('Thời gian (s)');
ylabel('Vị trí của xe (m)');
% title('Vị trí của xe con');
% legend('ADRC','ADRC+ZV','ADRC+ZVD','ADRC+ETM4');
legend('Mô phỏng','Thực nghiệm');
axis([0 10 0 0.25]);
xticks([0 2 4 6 8 10]);
% yticks([0 0.1 0.2]);
% legend('ADRC');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Dao dong cuoi thanh
figure('Position', B);
hold on 
grid on
H = xlsread('C:\Users\ADMIN\OneDrive - Hanoi University of Science and Technology\Tailieuhoctap\Stacker_crane\DATNCN\Simulation\DATA_THUCNGHIEM\adrc9.xlsx');
plot(H(:,1),H(:,7),'linewidth',1.2); 
xlabel('Thời gian (s)');
ylabel('Dao động (mm)');
legend('ADRC')
% title('Đồ thị dao động xét tại điểm cuối của thanh');
axis([0 10 -20 20]);
% xticks([0 2 4 6 8]);
% yticks([0 0.1 0.2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Phan tich pho
figure('Position', B);
hold on;
grid on;
N= xlsread('C:\Users\ADMIN\OneDrive - Hanoi University of Science and Technology\Tailieuhoctap\Stacker_crane\DATNCN\Simulation\DATA_THUCNGHIEM\ADRC\ADRC_dd.xlsx');
tmax=N(length(N(:,1)),1);
r=length(N(:,1));
dt=tmax/r;
Fs=1/dt;
time=0:dt:(r-2)*dt-dt;
l1=r;
fft_w=fft(N(:,3),l1)*(2/l1);
abs_w=abs(fft_w);
freq=0:(1/time(end)):Fs/2-(1/time(end));
plot(freq,abs_w(1:length(freq)),'LineWidth',1.2); 
xlabel('Frequency(Hz)');
ylabel('Amplitude(mm)');
% title('Phân tích phổ dao động xét tại điểm cuối của thanh');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%So sánh các kết quả chống rung điểm cuối
figure('Position', B);
%H = xlsread('C:\Users\Admin\Downloads\DATN\Du liệu\daodong.xlsx');
H = xlsread('C:\Users\ADMIN\OneDrive - Hanoi University of Science and Technology\Tailieuhoctap\Stacker_crane\DATNCN\Simulation\DATA_THUCNGHIEM\adrc9.xlsx');
plot(H(:,1),H(:,7),'b','linewidth',0.8);
hold on 
grid on 
plot(H(:,1),H(:,8),'r-.','linewidth',1);
plot(H(:,1),H(:,9),'k','linewidth',1);
plot(H(:,1),H(:,10),'m','linewidth',1.5);
xlabel('Thời gian (s)');
ylabel('Biên độ dao động (mm)');
legend('ADRC','ADRC+ZV','ADRC+ZVD','ADRC+ETM4');
title('So sánh kết quả các phương pháp chống rung');

figure('Position', [100 100 700 400]);
hold on 
grid on
subplot(2,2,1)
plot(H(:,1),H(:,7),'linewidth',1.2);
legend('ADRC');
xlabel('Thời gian (s)');
ylabel('ω(L,t)(m)')
axis([0 10 -20 20]);
grid on

subplot(2,2,2)
plot(H(:,1),H(:,8),'color',[0.12,0.66,0.12],'linewidth',1);
legend('ADRC+ZV');
xlabel('Thời gian (s)');
ylabel('ω(L,t)(m)')
axis([0 10 -20 20]);
grid on

subplot(2,2,3)
plot(H(:,1),H(:,9),'color',[0.85,0.33,0.10],'linewidth',1);
legend('ADRC+ZVD');
xlabel('Thời gian (s)');
ylabel('ω(L,t)(m)')
axis([0 10 -20 20]);
grid on

subplot(2,2,4)
plot(H(:,1),H(:,10),'k','linewidth',1.2);
legend('ADRC+ETM4');
xlabel('Thời gian (s)');
ylabel('ω(L,t)(m)')
axis([0 10 -20 20]);
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', B);
G = xlsread('C:\Users\ADMIN\OneDrive - Hanoi University of Science and Technology\Tailieuhoctap\Stacker_crane\DATNCN\Simulation\DATA_THUCNGHIEM\21052024\adrc1.xlsx');
hold on
tmax=G(length(G(:,1)),1);
r=length(G(:,1));
dt=tmax/r;
Fs=1/dt;
time=0:dt:(r-2)*dt-dt;
l1=r;
fft_w=fft(G(:,2),l1)*(2/l1);
abs_w=abs(fft_w);
freq=0:(1/time(end)):Fs/2-(1/time(end));
hold on
plot(freq,abs_w(1:length(freq)),'LineWidth',1.2); 
grid on
xlabel('Tần số (Hz)');
ylabel('Biên độ (mm)');

N= xlsread('C:\Users\ADMIN\OneDrive - Hanoi University of Science and Technology\Tailieuhoctap\Stacker_crane\DATNCN\Simulation\DATA_THUCNGHIEM\ADRC\ADRC_dd.xlsx');
tmax=N(length(N(:,1)),1);
r=length(N(:,1));
dt=tmax/r;
Fs=1/dt;
time=0:dt:(r-2)*dt-dt;
l1=r;
fft_w=fft(N(:,3),l1)*(2/l1);
abs_w=abs(fft_w);
freq=0:(1/time(end)):Fs/2-(1/time(end));
plot(freq,abs_w(1:length(freq)),'--','LineWidth',1.2); 
legend('Điểm giữa thang','Đỉnh thang')


% title('Phân tích phổ dao động xét tại điểm thứ 6 của thanh');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', B);
hold on 
grid on
K = xlsread('C:\Users\ADMIN\OneDrive - Hanoi University of Science and Technology\Tailieuhoctap\Stacker_crane\DATNCN\Simulation\DATA_THUCNGHIEM\adrc6.xlsx');
plot(K(:,1),K(:,7),'linewidth',1.2);
plot(H(:,1),H(:,7),'linewidth',1.2); 
xlabel('Thời gian (s)');
ylabel('Dao động (mm)');
legend('Điểm giữa thang','Đỉnh thang')
axis([0 10 -20 30]);
% title('Đồ thị dao động xét tại điểm thứ 6 của thanh');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', B);
% K = xlsread('C:\Users\Admin\Downloads\DATN\Du liệu\adrc6.xlsx');
plot(K(:,1),K(:,7),'b','linewidth',0.6);
hold on 
grid on 
plot(K(:,1),K(:,8),'r-.','linewidth',1);
plot(K(:,1),K(:,9),'k','linewidth',1);
plot(K(:,1),K(:,10),'m','linewidth',1.5);
xlabel('Thời gian (s)');
ylabel('Biên độ dao động (mm)');
legend('ADRC','ADRC+ZV','ADRC+ZVD','ADRC+ETM4');
title('So sánh kết quả các phương pháp chống rung');

figure('Position', [100 100 700 400]);
hold on 
grid on
subplot(2,2,1)
plot(K(:,1),K(:,7),'linewidth',1.2);
legend('ADRC');
xlabel('Thời gian (s)');
ylabel('ω(L/2,t)(m)')
axis([0 10 -10 10]);
grid on

subplot(2,2,2)
plot(K(:,1),K(:,8),'color',[0.12,0.66,0.12],'linewidth',1);
legend('ADRC+ZV');
xlabel('Thời gian (s)');
ylabel('ω(L/2,t)(m)')
axis([0 10 -10 10]);
grid on

subplot(2,2,3)
plot(K(:,1),K(:,9),'color',[0.85,0.33,0.10],'linewidth',1);
legend('ADRC+ZVD');
xlabel('Thời gian (s)');
ylabel('ω(L/2,t)(m)')
axis([0 10 -10 10]);
grid on

subplot(2,2,4)
plot(K(:,1),K(:,10),'k','linewidth',1.2);
legend('ADRC+ETM4');
xlabel('Thời gian (s)');
ylabel('ω(L/2,t)(m)')
axis([0 10 -10 10]);
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(8);
O= xlsread("C:\Users\ADMIN\OneDrive - Hanoi University of Science and Technology\Tailieuhoctap\Stacker_crane\DATNCN\Simulation\DATA_THUCNGHIEM\21052024\pid212.xlsx");
% hold on
% tmax=O(length(O(:,1)),1);
% r=length(O(:,1));
% dt=tmax/r;
% Fs=1/dt;
% time=0:dt:(r-2)*dt-dt;
% l1=r;
% fft_w=fft(O(:,2),l1)*(2/l1);
% abs_w=abs(fft_w);
% freq=0:(1/time(end)):Fs/2-(1/time(end));
% hold on
% plot(freq,abs_w(1:length(freq)),'LineWidth',0.8); 
% grid on
% xlabel('Tần số dao động (Hz)');
% ylabel('Biên độ tần số');
% title('Phân tích phổ dao động xét tại điểm cuối của thanh');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', B);
% O= xlsread('C:\Users\Admin\Downloads\DATN\Du liệu\21052024\pid212.xlsx');
plot(O(:,1),O(:,3),'b','linewidth',0.6);
hold on 
grid on 
legend('PID 212')
xlabel('Thời gian (s)');
ylabel('Biên độ dao động (mm)');
title('Đồ thị dao động xét tại điểm cuối của thanh');
