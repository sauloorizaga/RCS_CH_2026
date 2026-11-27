% =========================================================
% CPC_Table1.m
% Reproduces Table 1 and Figure 3 of the paper:
%   L2 errors and convergence rates at T_f = 1 for
%   CS (1st order), BDF2-CS (2nd order) and RCS (RCSX-iter2).
% GPU version (MATLAB gpuArray).
%
% Needs:  U_trueIC.mat  (smooth initial condition, variable U_ic)
%         U_ref.mat     (reference solution, variable Uexact), or set
%                       RUN_REFERENCE = true to build it (slow).
%
% S. Orizaga, "A Second-Order Richardson--Convex Splitting Method for the
% Cahn--Hilliard Equation: Stability Analysis and GPU-Accelerated 3D
% Computations", Computer Physics Communications.
% =========================================================
clc; clear;

% =========================
% PARAMETERS (paper values)
% =========================
M = 1; a = -M*pi; b = M*pi;
N = 256;
epsilon = 0.1; eps2 = epsilon^2;
a_split = 2;            % splitting parameter used in Table 1
tfinal  = 1.0;
m_steps = 5;            % CS substeps used to compute U^1 for RCS
RUN_REFERENCE = false;  % true: build U_ref.mat (dt_ref = 1e-5), false: load it

dt_values = [0.1 0.04 0.02 0.01 0.005 0.0025 0.0025/2 0.0025/4 0.0025/8];

% =========================
% GRID AND WAVENUMBERS
% =========================
h = (b-a)/N;
x = [a+h:h:b];
[X,Y] = meshgrid(x,x);
k = [[0:N/2] [-N/2+1:-1]]./(M);
[kx,ky] = meshgrid(k.^2,k.^2);
k2 = kx + ky;
k4 = k2.^2;
k2g = gpuArray(k2);
k4g = gpuArray(k4);

% =========================
% INITIAL CONDITION
% =========================
S = load('U_trueIC.mat');
Uinit = S.U_ic;

% =========================
% REFERENCE SOLUTION (RCSX-iter2, dt_ref = 1e-5)
% =========================
if RUN_REFERENCE
    dt_ref = 1e-5;
    lhs_ref      = gpuArray(1 + dt_ref*(eps2*k4 + a_split*k2));
    lhs_ref_fine = gpuArray(1 + (dt_ref/2)*(eps2*k4 + a_split*k2));
    Uref = gpuArray(Uinit);
    hat_U = fft2(Uref);
    % 2-step bootstrap for U1
    m_boot = 2; dt_boot = dt_ref/m_boot;
    lhs_boot = gpuArray(1 + dt_boot*(eps2*k4 + a_split*k2));
    hat_U_prev = hat_U; U_tmp = Uref;
    for m = 1:m_boot
        fU_tmp = U_tmp.^3 - (1+a_split)*U_tmp;
        U_tmp = real(ifft2((fft2(U_tmp) + dt_boot*(-k2g.*fft2(fU_tmp))) ./ lhs_boot));
    end
    hat_U = fft2(U_tmp);
    t = dt_ref;
    while t < tfinal - dt_ref/2
        hat_U_n = hat_U;
        U_n   = real(ifft2(hat_U_n));
        U_nm1 = real(ifft2(hat_U_prev));
        U_extrap = 2*U_n - U_nm1;
        U_c = U_extrap;
        for i = 1:2
            U_c = real(ifft2((hat_U_n + dt_ref*(-k2g.*fft2(U_c.^3 - (1+a_split)*U_c))) ./ lhs_ref));
        end
        U_h = U_extrap;
        for i = 1:2
            U_h = real(ifft2((hat_U_n + (dt_ref/2)*(-k2g.*fft2(U_h.^3 - (1+a_split)*U_h))) ./ lhs_ref_fine));
        end
        U_f = U_h; hat_Uh = fft2(U_h);
        for i = 1:2
            U_f = real(ifft2((hat_Uh + (dt_ref/2)*(-k2g.*fft2(U_f.^3 - (1+a_split)*U_f))) ./ lhs_ref_fine));
        end
        hat_U_prev = hat_U_n;
        hat_U = fft2(2*U_f - U_c);
        t = t + dt_ref;
    end
    Uexact = gather(real(ifft2(hat_U)));
    save('U_ref.mat','Uexact');
else
    load('U_ref.mat','Uexact');
end
fprintf('Reference solution ready.\n');

% =========================
% MAIN LOOP
% =========================
err_CS   = zeros(size(dt_values));
err_BDF2 = zeros(size(dt_values));
err_RCS  = zeros(size(dt_values));

for idx = 1:length(dt_values)
    dt = dt_values(idx);
    lhs   = gpuArray(1 + dt*(eps2*k4 + a_split*k2));
    lhs_f = gpuArray(1 + (dt/2)*(eps2*k4 + a_split*k2));

    % ---------- CS (1st order) ----------
    U = gpuArray(Uinit); hat_U = fft2(U); t = 0;
    while t < tfinal - dt/2
        fU = U.^3 - (1+a_split)*U;
        hat_U = (hat_U + dt*(-k2g.*fft2(fU))) ./ lhs;
        U = real(ifft2(hat_U));
        t = t + dt;
    end
    err_CS(idx) = sqrt(h^2*sum(sum((gather(U) - Uexact).^2)));

    % ---------- BDF2-CS (classical, extrapolated nonlinearity) ----------
    lhs_bdf2 = gpuArray(1.5 + dt*(eps2*k4 + a_split*k2));
    U0 = gpuArray(Uinit);
    hat_bdf2_nm1 = fft2(U0);
    hat_bdf2_n = (hat_bdf2_nm1 + dt*(-k2g.*fft2(U0.^3 - (1+a_split)*U0))) ./ lhs;   % one CS step for U^1
    t = dt;
    while t < tfinal - dt/2
        U_extrap = 2*real(ifft2(hat_bdf2_n)) - real(ifft2(hat_bdf2_nm1));
        rhs = 2*hat_bdf2_n - 0.5*hat_bdf2_nm1 + dt*(-k2g.*fft2(U_extrap.^3 - (1+a_split)*U_extrap));
        hat_bdf2_nm1 = hat_bdf2_n;
        hat_bdf2_n = rhs ./ lhs_bdf2;
        t = t + dt;
    end
    err_BDF2(idx) = sqrt(h^2*sum(sum((gather(real(ifft2(hat_bdf2_n))) - Uexact).^2)));

    % ---------- RCS (RCSX-iter2: predictor 2U^n - U^{n-1}, two sweeps per substep) ----------
    U = gpuArray(Uinit); hat_U = fft2(U);
    hat_U_prev = hat_U;                                   % U^0
    dt_tiny = dt/m_steps;                                 % U^1 from m_steps CS substeps
    lhs_tiny = gpuArray(1 + dt_tiny*(eps2*k4 + a_split*k2));
    for m = 1:m_steps
        U = real(ifft2((fft2(U) + dt_tiny*(-k2g.*fft2(U.^3 - (1+a_split)*U))) ./ lhs_tiny));
    end
    hat_U = fft2(U); t = dt;
    while t < tfinal - dt/2
        hat_U_n = hat_U;
        U_extrap = 2*U - real(ifft2(hat_U_prev));
        % Coarse step (dt)
        Uc = U_extrap;
        for i = 1:2
            Uc = real(ifft2((hat_U_n + dt*(-k2g.*fft2(Uc.^3 - (1+a_split)*Uc))) ./ lhs));
        end
        % Fine step 1 (dt/2)
        Uh = U_extrap;
        for i = 1:2
            Uh = real(ifft2((hat_U_n + (dt/2)*(-k2g.*fft2(Uh.^3 - (1+a_split)*Uh))) ./ lhs_f));
        end
        % Fine step 2 (dt/2)
        Uf = Uh; hat_Uh = fft2(Uh);
        for i = 1:2
            Uf = real(ifft2((hat_Uh + (dt/2)*(-k2g.*fft2(Uf.^3 - (1+a_split)*Uf))) ./ lhs_f));
        end
        % Richardson combination
        U = 2*Uf - Uc;
        hat_U_prev = hat_U_n;
        hat_U = fft2(U);
        t = t + dt;
    end
    err_RCS(idx) = sqrt(h^2*sum(sum((gather(U) - Uexact).^2)));

    fprintf('dt=%-10g CS=%.3e  BDF2-CS=%.3e  RCS=%.3e\n', dt, err_CS(idx), err_BDF2(idx), err_RCS(idx));
end

% =========================
% TABLE 1 (rows dt = 0.04 ... 0.0003125)
% =========================
ord = @(e) [NaN, log2(e(1:end-1)./e(2:end))];
oCS = ord(err_CS); oB = ord(err_BDF2); oR = ord(err_RCS);
fprintf('\n%-10s | %-10s %-6s | %-10s %-6s | %-10s %-6s\n', 'dt', 'CS', 'Order', 'BDF2-CS', 'Order', 'RCS', 'Order');
fprintf('%s\n', repmat('-', 1, 72));
for i = 2:length(dt_values)
    if i == 2
        fprintf('%-10g | %-10.2e %-6s | %-10.2e %-6s | %-10.2e %-6s\n', dt_values(i), ...
            err_CS(i), '---', err_BDF2(i), '---', err_RCS(i), '---');
    else
        fprintf('%-10g | %-10.2e %-6.2f | %-10.2e %-6.2f | %-10.2e %-6.2f\n', dt_values(i), ...
            err_CS(i), oCS(i), err_BDF2(i), oB(i), err_RCS(i), oR(i));
    end
end
save('CPC_Table1_errors.mat', 'dt_values', 'err_CS', 'err_BDF2', 'err_RCS');

% =========================
% FIGURE 3: error plot
% =========================
figure(12); clf;
set(gcf, 'Color', 'w', 'Position', [100, 100, 800, 600]);
loglog(dt_values, err_CS, 'o-', 'Color', [0.6 0.6 0.6], ...
    'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Standard CS');
hold on;
loglog(dt_values, err_BDF2, 'ks-', 'LineWidth', 2, ...
    'MarkerSize', 8, 'DisplayName', 'BDF2--CS');
loglog(dt_values, err_RCS, 'o-', 'Color', [0 0.45 0.74], ...
    'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', [0 0.45 0.74], ...
    'DisplayName', 'RCS');
dt_ref = dt_values(end-4:end);
loglog(dt_ref, .5 * err_CS(end) * (dt_ref/dt_ref(end)), ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(dt_ref(2), .6 * err_CS(end-2), '$O(\Delta t)$', 'Interpreter', 'latex', 'FontSize', 14, 'Color', [0.4 0.4 0.4]);
loglog(dt_ref, 70 * err_RCS(end) * (dt_ref/dt_ref(end)).^2, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(dt_ref(2), 165 * err_RCS(end-2), '$O(\Delta t^2)$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
grid on; grid minor;
ax = gca; ax.FontSize = 16;
ax.TickLabelInterpreter = 'latex';
xlabel('Time Step $\Delta t$', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('$L^2$ Error', 'FontSize', 20, 'Interpreter', 'latex');
legend('Location', 'southeast', 'FontSize', 15, 'Interpreter', 'latex', 'Box', 'off');
axis tight;
