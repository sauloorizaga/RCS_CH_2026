% =========================================================
% Point 5: BDF2-CS vs RCSX-iter2 at T = 100 for several dt
% Same IC, eps, grid; final morphologies (2x4) + domain table
% =========================================================
clc; clear; close all;

N = 256; tfinal = 100; epsilon = 0.05; eps2 = epsilon^2;
a_split = 2; m_boot = 5;
dts = [0.01 0.005 0.0025 0.00125];

M = 1; a = -M*pi; b = M*pi; h = (b-a)/N;
x = gpuArray(linspace(a+h, b, N)); [X, Y] = meshgrid(x, x);
k = gpuArray([[0:N/2] [-N/2+1:-1]] ./ M);
[k1x, k1y] = meshgrid(k, k); [kx, ky] = meshgrid(k.^2, k.^2);
k2g = kx + ky; k4g = k2g.^2;

rng(100,'twister'); Uinit = gpuArray(0.05*rand(N, N) - 0.025);

Energy = @(U) gather(h^2*sum(sum( eps2/2*(real(ifft2(1i*k1x.*fft2(U))).^2 + ...
         real(ifft2(1i*k1y.*fft2(U))).^2) + 0.25*(U.^2-1).^2 )));

nd = numel(dts);
UB = cell(nd,1); UR = cell(nd,1);
EB = zeros(nd,1); ER = EB; DB = EB; DR = EB;

for j = 1:nd
    dt = dts(j);
    lhs_cs = 1 + dt*(eps2*k4g + a_split*k2g);
    lhs_f  = 1 + (dt/2)*(eps2*k4g + a_split*k2g);
    lhs_bdf2 = 1.5 + dt*(eps2*k4g + a_split*k2g);

    % BDF2 start: one CS step
    hb_nm1 = fft2(Uinit);
    hb_n = (hb_nm1 + dt*(-k2g.*fft2(Uinit.^3-(1+a_split)*Uinit))) ./ lhs_cs;

    % RCS start: m_boot tiny CS steps
    U_boot = Uinit; dt_tiny = dt/m_boot;
    lhs_tiny = 1 + dt_tiny*(eps2*k4g + a_split*k2g);
    for m = 1:m_boot
        U_boot = real(ifft2((fft2(U_boot) + dt_tiny*(-k2g.*fft2(U_boot.^3-(1+a_split)*U_boot))) ./ lhs_tiny));
    end
    hr_n = fft2(U_boot); hr_nm1 = fft2(Uinit);

    t = dt;
    while t < tfinal + dt/2
        % RCSX-iter2
        Un = real(ifft2(hr_n)); Ue = 2*Un - real(ifft2(hr_nm1));
        Uc = Ue; for i=1:2, Uc = real(ifft2((hr_n + dt*(-k2g.*fft2(Uc.^3-(1+a_split)*Uc))) ./ lhs_cs)); end
        Uh = Ue; for i=1:2, Uh = real(ifft2((hr_n + (dt/2)*(-k2g.*fft2(Uh.^3-(1+a_split)*Uh))) ./ lhs_f)); end
        Uf = Uh; hUh = fft2(Uh);
        for i=1:2, Uf = real(ifft2((hUh + (dt/2)*(-k2g.*fft2(Uf.^3-(1+a_split)*Uf))) ./ lhs_f)); end
        U_r = 2*Uf - Uc;
        % BDF2-CS classical
        Ue_b = 2*real(ifft2(hb_n)) - real(ifft2(hb_nm1));
        hb_np1 = (2*hb_n - 0.5*hb_nm1 + dt*(-k2g.*fft2(Ue_b.^3-(1+a_split)*Ue_b))) ./ lhs_bdf2;
        % update
        hr_nm1 = hr_n; hr_n = fft2(U_r);
        hb_nm1 = hb_n; hb_n = hb_np1;
        t = t + dt;
    end

    UR{j} = gather(real(ifft2(hr_n))); UB{j} = gather(real(ifft2(hb_n)));
    ER(j) = Energy(gpuArray(UR{j}));   EB(j) = Energy(gpuArray(UB{j}));
    DR(j) = count_periodic(UR{j} > 0); DB(j) = count_periodic(UB{j} > 0);
    fprintf('dt=%-8g BDF2: %d domains, E=%.4f | RCS: %d domains, E=%.4f\n', ...
        dt, DB(j), EB(j), DR(j), ER(j));
end
save('point5_morphology.mat','dts','UB','UR','EB','ER','DB','DR');

% ---------- Table ----------
fprintf('\n%-10s | %-8s %-9s | %-8s %-9s\n','dt','Dom BDF2','E BDF2','Dom RCS','E RCS');
for j = 1:nd
    fprintf('%-10g | %-8d %-9.4f | %-8d %-9.4f\n', dts(j), DB(j), EB(j), DR(j), ER(j));
end

% ---------- 2x4 figure ----------
Xg = gather(X); Yg = gather(Y);
figure(40); clf; set(gcf,'Color','w','Position',[50 50 1500 700]);
for j = 1:nd
    subplot(2,nd,j);
    pcolor(Xg, Yg, UB{j}); shading interp; axis equal tight off;
    title(sprintf('BDF2-CS, \\Delta t = %g\n%d domains, E = %.3f', dts(j), DB(j), EB(j)), 'FontSize', 11);
    subplot(2,nd,nd+j);
    pcolor(Xg, Yg, UR{j}); shading interp; axis equal tight off;
    title(sprintf('RCS, \\Delta t = %g\n%d domains, E = %.3f', dts(j), DR(j), ER(j)), 'FontSize', 11);
end
exportgraphics(gcf, 'MorphologyRefinement.png', 'Resolution', 300);

% ================= function =================
function n = count_periodic(mask)
% number of connected regions of mask on a periodic grid (4-connectivity),
% without toolboxes (flood fill with periodic neighbors)
mask = logical(mask); [Nr, Nc] = size(mask);
seen = false(Nr, Nc); n = 0;
queue = zeros(Nr*Nc, 2);
for i0 = 1:Nr
    for j0 = 1:Nc
        if mask(i0,j0) && ~seen(i0,j0)
            n = n + 1;
            head = 1; tail = 1; queue(1,:) = [i0 j0]; seen(i0,j0) = true;
            while head <= tail
                i = queue(head,1); j = queue(head,2); head = head + 1;
                nb = [mod(i-2,Nr)+1 j; mod(i,Nr)+1 j; i mod(j-2,Nc)+1; i mod(j,Nc)+1];
                for q = 1:4
                    ii = nb(q,1); jj = nb(q,2);
                    if mask(ii,jj) && ~seen(ii,jj)
                        seen(ii,jj) = true; tail = tail + 1; queue(tail,:) = [ii jj];
                    end
                end
            end
        end
    end
end
end