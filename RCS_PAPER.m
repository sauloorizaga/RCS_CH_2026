% =========================================================
% BDF2-CS vs RCSX-iter2: 2D Long-Run Challenge (Side-by-Side)
% =========================================================
clc; clear; close all;

% --- Configuration ---
N = 256;             % Resolution
tfinal = 100.0;        % Target for phase separation dynamics
dt = .01;           % Time step
epsilon = 0.05;      % Interfacial width
eps2 = epsilon^2;
a_split = 2;       % Splitting parameter (Theoretical a >= 1 for RCS)
m_boot = 50;          % Tiny steps for bootstrap

% --- Domain and Operators ---
M = 1; a = -M*pi; b = M*pi; h = (b-a)/N;
x = gpuArray(linspace(a+h, b, N));
[X, Y] = meshgrid(x, x);
k = gpuArray([[0:N/2] [-N/2+1:-1]] ./ M);
[k1x, k1y] = meshgrid(k, k);
[kx, ky] = meshgrid(k.^2, k.^2);
k2g = kx + ky; k4g = k2g.^2;

% --- Initial Condition ---
rng(100,'twister')
%rng(100,'twister'); % For reproducibility
Uinit = 0.05*rand(N, N) - 0.025; % Clean random start
Uinit=gpuArray(Uinit)

% --- LHS Operators ---
lhs_cs = 1 + dt*(eps2*k4g + a_split*k2g);
lhs_f = 1 + (dt/2)*(eps2*k4g + a_split*k2g);
lhs_bdf2 = 1.5 + dt*(eps2*k4g + a_split*k2g);

% --- Method 1: BDF2 State Initialization ---
hat_bdf2_nm1 = fft2(Uinit);
% Initial CS step for BDF2 bootstrap
fU_init = Uinit.^3 - (1+a_split)*Uinit;
hat_bdf2_n = (hat_bdf2_nm1 + dt*(-k2g.*fft2(fU_init))) ./ lhs_cs;

% --- Method 2: RCSX State Initialization ---
% Improved Bootstrap (tiny steps)
U_boot = Uinit;
dt_tiny = dt / m_boot;
lhs_tiny = 1 + dt_tiny*(eps2*k4g + a_split*k2g);
for m = 1:m_boot
    fU_boot = U_boot.^3 - (1+a_split)*U_boot;
    hat_U_boot = (fft2(U_boot) + dt_tiny*(-k2g.*fft2(fU_boot))) ./ lhs_tiny;
    U_boot = real(ifft2(hat_U_boot));
end
hat_rcsx_n = fft2(U_boot);
hat_rcsx_nm1 = fft2(Uinit); % For initial extrapolation

% --- Pre-allocate Energy and Time ---
steps = ceil(tfinal/dt);
E_bdf2 = zeros(steps, 1); E_rcsx = zeros(steps, 1);
time_axis = zeros(steps, 1);

% --- Simulation Loop ---
t = dt; it = 1;
figure(1); set(gcf, 'Position', [100, 100, 1200, 500]);

while t < tfinal + dt/2
    % ---------------------------------
    % RCSX-iter2 Logic
    % ---------------------------------
    Un_rc = real(ifft2(hat_rcsx_n));
    Unm1_rc = real(ifft2(hat_rcsx_nm1));
    U_extrap = 2*Un_rc - Unm1_rc; % High-order guess
    
    % Coarse solve (2 iters)
    Uc = U_extrap;
    for i = 1:2
        fUc = Uc.^3 - (1+a_split)*Uc;
        Uc = real(ifft2((hat_rcsx_n + dt*(-k2g.*fft2(fUc))) ./ lhs_cs));
    end
    % Fine solve (2 iters per half-step)
    Uh = U_extrap;
    for i = 1:2
        fUh = Uh.^3 - (1+a_split)*Uh;
        Uh = real(ifft2((hat_rcsx_n + (dt/2)*(-k2g.*fft2(fUh))) ./ lhs_f));
    end
    Uf = Uh; hat_Uh = fft2(Uh);
    for i = 1:2
        fUf = Uf.^3 - (1+a_split)*Uf;
        Uf = real(ifft2((hat_Uh + (dt/2)*(-k2g.*fft2(fUf))) ./ lhs_f));
    end
    U_rcsx_new = 2*Uf - Uc; % Richardson Mic Drop
    
    % ---------------------------------
    % BDF2-CS Logic
    % ---------------------------------
    Un_b = real(ifft2(hat_bdf2_n));
    Unm1_b = real(ifft2(hat_bdf2_nm1));
    U_ext_b = 2*Un_b - Unm1_b;
    fhat_bdf2 = fft2(U_ext_b.^3 - (1+a_split)*U_ext_b);
    rhs_bdf2 = 2*hat_bdf2_n - 0.5*hat_bdf2_nm1 + dt*(-k2g.*fhat_bdf2);
    hat_bdf2_np1 = rhs_bdf2 ./ lhs_bdf2;
    
    % ---------------------------------
    % Energy & Diagnostics
    % ---------------------------------
    % Precise Energy = Grad term (Spectral) + Potential
    E_rcsx(it) = gather(h^2 * sum(sum( eps2/2 * (real(ifft2(1i*k1x.*fft2(U_rcsx_new))).^2 + ...
                 real(ifft2(1i*k1y.*fft2(U_rcsx_new))).^2) + 0.25*(U_rcsx_new.^2 - 1).^2 )));
    
    U_bdf2_plot = real(ifft2(hat_bdf2_np1));
    E_bdf2(it) = gather(h^2 * sum(sum( eps2/2 * (real(ifft2(1i*k1x.*hat_bdf2_np1)).^2 + ...
                 real(ifft2(1i*k1y.*hat_bdf2_np1)).^2) + 0.25*(U_bdf2_plot.^2 - 1).^2 )));
    
    time_axis(it) = t;

    % Visual Evolution (Every 10 steps to keep it smooth)
    if mod(it, 500) == 0
        subplot(1,2,1); pcolor(X, Y, gather(U_bdf2_plot)); shading interp; axis equal tight off;
        title(['BDF2-CS | t = ' num2str(t, '%.2f')]);
        
        subplot(1,2,2); pcolor(X, Y, gather(U_rcsx_new)); shading interp; axis equal tight off;
        title(['RCSX-iter2 (High Fidelity) | t = ' num2str(t, '%.2f')]);
        drawnow;
    end
    
    % Update states
    hat_rcsx_nm1 = hat_rcsx_n; hat_rcsx_n = fft2(U_rcsx_new);
    hat_bdf2_nm1 = hat_bdf2_n; hat_bdf2_n = hat_bdf2_np1;
    t = t + dt; it = it + 1;
end

% --- Final Energy Comparison ---
figure(20);
% subplot(1,2,1);
% plot(time_axis(1:it-1), E_bdf2(1:it-1), 'r', 'LineWidth', 2); hold on;
% plot(time_axis(1:it-1), E_rcsx(1:it-1), 'b--', 'LineWidth', 2);
% xlabel('Time'); ylabel('Free Energy'); title('Energy Dissipation');
% legend('BDF2-CS', 'RCSX-iter2'); grid on;

% subplot(1,2,2);
semilogx(time_axis(1:it-1), E_bdf2(1:it-1), 'r', 'LineWidth', 2); hold on;
semilogx(time_axis(1:it-1), E_rcsx(1:it-1), 'b--', 'LineWidth', 2);
xlabel('Time (log)'); ylabel('Free Energy'); title('Energy (Log-Scale)');
grid on;legend('BDF2-CS', 'RCSX-iter2');
