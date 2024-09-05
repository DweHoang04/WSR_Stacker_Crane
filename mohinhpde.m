%clc; clear; close all
function stacker
    % Các tham số mô hình
    pA = 2.4882; % mật độ của thanh kg/m
    EI = 1.05; % độ cứng uốn (flexural rigidity) Nm^2
    mw = 13.1; % khối lượng xe con
    mk = 0.5; % khối lượng đỉnh thang
    L = 0.6; % Chiều dài của thanh
    F = 10; % Lực tác dụng
    N = 0;
    n = 20; % Số lượng phần tử
    m = 0;
    
    x = linspace(0, L, n-1);
    t = linspace(0, 10, 10);

    % Sử dụng hàm vô danh để truyền tham số vào các hàm yêu cầu bởi pdepe
    pde = @(x,t,u,DuDx) beampde(x,t,u,DuDx, EI, pA);
    bc = @(xl,ul,xr,ur,t) beambc(xl,ul,xr,ur,t, EI, mw, mk, F, L);
    ic = @(x) beamic(x, L);

    sol = pdepe(m, pde, ic, bc, x, t);

    % Vẽ kết quả
    w = sol(:,:,1);
    figure; plot(t, w); grid;
    xlabel('Time'); ylabel('Displacement'); 
    title('Displacement at Beam');

    figure; plot(t,sol(:,n-1,1)); grid;
    xlabel('Time'); ylabel('Displacement'); 
    title('Beam Displacement at the final point');

    figure; plot(t,sol(:,1,1)); grid;
    xlabel('Time'); ylabel('Displacement'); 
    title('Beam Displacement at the first point');
end

function [c,f,s] = beampde(x,t,u,DuDx, EI, pA)
    % pA * đạo hàm bậc hai theo thời gian của u
    c = pA * diff(u, t, 2); % Nhân pA với đạo hàm bậc hai của u theo thời gian
    f = EI * DuDx(3);  % Đây là EI * đạo hàm bậc ba của u theo x
    s = 0;  % Điều này thể hiện không có phần nguồn hoặc thành phần tự do
end

function u0 = beamic(x, L)
    % Điều kiện ban đầu
    u0 = 0;
end
	
function [pl,ql,pr,qr] = beambc(xl,ul,xr,ur,t, EI, mw, mk, F, L)
    % Điều kiện biên
    pl = mw * diff(ul, t, 2) - F + EI * subs(diff(ul, xl, 3), xl, 0);
    ql = 0;
    pr = mk * diff(ur, t, 2) - EI * subs(diff(ur, xr, 3), xr, L);
    qr = 0;
end
