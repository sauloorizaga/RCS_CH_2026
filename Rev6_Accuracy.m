% ===== Point 6 (2D): error and runtime, RCSX-iter2 vs BDF2-CS classical =====
clc; clear;
N = 256; M = 1; a = -M*pi; b = M*pi; h = (b-a)/N;
k = [[0:N/2] [-N/2+1:-1]]./M; [kx,ky] = meshgrid(k.^2,k.^2);
k2g = gpuArray(kx+ky); k4g = k2g.^2;
epsilon = 0.1; eps2 = epsilon^2; a_split = 2; tfinal = 1.0; m_boot = 50;
S = load('U_trueIC.mat'); Uinit = gpuArray(S.U_ic);
g = gpuDevice;

% ---------- Reference: RCSX-iter2 at dt = 1e-5 ----------
dt_ref = 1e-5;
Uexact = gather(rcsx_iter2(Uinit, dt_ref, tfinal, k2g, k4g, eps2, a_split, m_boot));
fprintf('Reference built.\n');

% ---------- Warm-up (GPU compilation), not timed ----------
rcsx_iter2(Uinit, 0.04, 0.08, k2g, k4g, eps2, a_split, m_boot);
bdf2_classic(Uinit, 0.04, 0.08, k2g, k4g, eps2, a_split);
wait(g);

dts = [0.04 0.02 0.01 0.005 0.0025 0.00125 0.000625 0.0003125];
errR = zeros(size(dts)); errB = errR; timeR = errR; timeB = errR;
for j = 1:numel(dts)
    dt = dts(j);
    wait(g); tic;
    U = rcsx_iter2(Uinit, dt, tfinal, k2g, k4g, eps2, a_split, m_boot);
    wait(g); timeR(j) = toc;
    errR(j) = sqrt(h^2*sum(sum((gather(U)-Uexact).^2)));

    wait(g); tic;
    U = bdf2_classic(Uinit, dt, tfinal, k2g, k4g, eps2, a_split);
    wait(g); timeB(j) = toc;
    errB(j) = sqrt(h^2*sum(sum((gather(U)-Uexact).^2)));

    fprintf('dt=%-9g RCS: err=%.3e time=%.3fs | BDF2: err=%.3e time=%.3fs\n', ...
        dt, errR(j), timeR(j), errB(j), timeB(j));
end
save('point6_2D_cost.mat','dts','errR','errB','timeR','timeB');

% ---------- Table at equal target error (log-log interpolation) ----------
targets = [1e-2 1e-3 1e-4];
fprintf('\n%-10s | %-10s %-8s %-9s | %-10s %-8s %-9s\n', ...
    'target','dt RCS','steps','time RCS','dt BDF2','steps','time BDF2');
for e = targets
    dR = 10^interp1(log10(errR), log10(dts),   log10(e));
    tR = 10^interp1(log10(errR), log10(timeR), log10(e));
    dB = 10^interp1(log10(errB), log10(dts),   log10(e));
    tB = 10^interp1(log10(errB), log10(timeB), log10(e));
    fprintf('%-10.0e | %-10.2e %-8d %-9.3f | %-10.2e %-8d %-9.3f\n', ...
        e, dR, round(tfinal/dR), tR, dB, round(tfinal/dB), tB);
end

% ---------- Plot: error versus runtime ----------
figure(30); clf; set(gcf,'Color','w');
loglog(timeB, errB, 'ks-', 'LineWidth', 2, 'MarkerSize', 8); hold on;
loglog(timeR, errR, 'bp-', 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', 'b');
grid on; ax = gca; ax.FontSize = 14;
xlabel('Runtime (s)', 'FontSize', 16); ylabel('L^2 error', 'FontSize', 16);
legend('BDF2-CS', 'RCSX-iter2', 'Location', 'NorthEast');
title('Error versus runtime (N = 256, T_f = 1)', 'FontSize', 16);

% ================= functions =================
function U = rcsx_iter2(U0, dt, tfinal, k2g, k4g, eps2, a, m_boot)
    lhs = 1+dt*(eps2*k4g+a*k2g); lhs_f = 1+(dt/2)*(eps2*k4g+a*k2g);
    dtt = dt/m_boot; lt = 1+dtt*(eps2*k4g+a*k2g); Ub = U0;
    for m = 1:m_boot
        Ub = real(ifft2((fft2(Ub)+dtt*(-k2g.*fft2(Ub.^3-(1+a)*Ub)))./lt));
    end
    Unm1 = U0; U = Ub; nsteps = round(tfinal/dt);
    for it = 2:nsteps
        hUn = fft2(U); Ue = 2*U-Unm1;
        Uc = Ue; for i=1:2, Uc = real(ifft2((hUn+dt*(-k2g.*fft2(Uc.^3-(1+a)*Uc)))./lhs)); end
        Uh = Ue; for i=1:2, Uh = real(ifft2((hUn+(dt/2)*(-k2g.*fft2(Uh.^3-(1+a)*Uh)))./lhs_f)); end
        Uf = Uh; hUh = fft2(Uh);
        for i=1:2, Uf = real(ifft2((hUh+(dt/2)*(-k2g.*fft2(Uf.^3-(1+a)*Uf)))./lhs_f)); end
        Unm1 = U; U = 2*Uf-Uc;
    end
end

function U = bdf2_classic(U0, dt, tfinal, k2g, k4g, eps2, a)
    lhs = 1+dt*(eps2*k4g+a*k2g); lhs_b = 1.5+dt*(eps2*k4g+a*k2g);
    h_nm1 = fft2(U0);
    h_n = (h_nm1+dt*(-k2g.*fft2(U0.^3-(1+a)*U0)))./lhs;
    nsteps = round(tfinal/dt);
    for it = 2:nsteps
        Ue = 2*real(ifft2(h_n))-real(ifft2(h_nm1));
        h_np1 = (2*h_n-0.5*h_nm1+dt*(-k2g.*fft2(Ue.^3-(1+a)*Ue)))./lhs_b;
        h_nm1 = h_n; h_n = h_np1;
    end
    U = real(ifft2(h_n));
end