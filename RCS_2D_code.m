% RCS 2D Solver for the Cahn--Hilliard Equation
% Self-contained GPU implementation
% Based on the RCSX-iter2 3D solver
clc
clear
close all
tic;
% --- Physical and numerical parameters ---
a = 0;b = 2*pi;
N = 256;
Tf = 50;
dt = 0.01;
eps = 0.1;
eps2 = eps^2;
a_split = 1;
h = (b-a)/N;
% --- GPU Mesh and Spectral Operators ---
x = gpuArray(linspace(a, b-h, N));
k = gpuArray([[0:N/2] [-N/2+1:-1]]);
[Kx, Ky] = meshgrid(k, k);
k2g = Kx.^2 + Ky.^2;k4g = k2g.^2;
% --- Linear operators ---
lhs_cs = 1 + dt*(eps2*k4g + a_split*k2g);
lhs_f  = 1 + (dt/2)*(eps2*k4g + a_split*k2g);
% --- Initial condition ---
rng(1527, 'twister');
U = gpuArray(0.05*rand(N,N) - 0.025);
U_old = U;hat_U_n = fft2(U);
InitialMass = gather(sum(U(:)));
% --- Time setup ---
steps = ceil(Tf/dt);
E_history = zeros(steps,1);
time_axis = zeros(steps,1);
fprintf('Starting 2D RCSX-iter2 (N=%d, dt=%.3f)...\n', N, dt);
fprintf('GPU computation in progress...\n');
tic;t = 0;
it = 1;
while t < Tf - dt*0.01
    % -------------------------------------------------------------
    % Richardson extrapolation
    % -------------------------------------------------------------
    U_extrap = 2*U - U_old;
    % --- Coarse step ---
    Uc = U_extrap;
    for i = 1:2
        fUc = Uc.^3 - (1+a_split)*Uc;
        Uc = real(ifft2( ...
            (hat_U_n + dt*(-k2g.*fft2(fUc))) ./ lhs_cs));
    end
    % --- Half step ---
    Uh = U_extrap;
    for i = 1:2
        fUh = Uh.^3 - (1+a_split)*Uh;
        Uh = real(ifft2( ...
            (hat_U_n + (dt/2)*(-k2g.*fft2(fUh))) ./ lhs_f));
    end
    % --- Second half step ---
    Uf = Uh;
    hat_Uh = fft2(Uh);
    for i = 1:2
        fUf = Uf.^3 - (1+a_split)*Uf;
        Uf = real(ifft2( ...
            (hat_Uh + (dt/2)*(-k2g.*fft2(fUf))) ./ lhs_f));
    end
    % --- Richardson extrapolation ---
    U_next = 2*Uf - Uc;
    % -------------------------------------------------------------
    % Free energy
    % -------------------------------------------------------------
    hatU_next = fft2(U_next);
    E_grad = (eps2/2) * h^2 * ...
        sum(k2g(:).*abs(hatU_next(:)).^2) / N^2;
    E_pot = 0.25 * h^2 * ...
        sum((U_next(:).^2 - 1).^2);
    E_history(it) = gather(E_grad + E_pot);
    % -------------------------------------------------------------
    % Update
    % -------------------------------------------------------------
    U_old = U;
    U = U_next;
    hat_U_n = fft2(U);
    t = t + dt;time_axis(it) = t;
    it = it + 1;
end
computing_time = toc;
% -------------------------------------------------------------
% Final solution and mass conservation
% -------------------------------------------------------------
U = gather(U);
FinalMass = sum(U(:));
RealMassError = abs(InitialMass - FinalMass) * h^2;
fprintf('\nSimulation Completed.\n');
fprintf('Final Mass: %.16e\n', FinalMass);
fprintf('Mass Error Computed: %.16e\n', RealMassError);
minutes = computing_time/60;
hours = computing_time/3600;
fprintf('Computing Time: %.4f seconds\n', computing_time);
fprintf('Computing Time: %.4f minutes\n', minutes);
fprintf('Computing Time: %.4f hours\n', hours);
% -------------------------------------------------------------
% 2D Solution
% -------------------------------------------------------------
figure(1); clf;
set(gcf, 'Color', 'w');
x_cpu = linspace(a,b,N);
[X,Y] = meshgrid(x_cpu,x_cpu);
pcolor(X,Y,U);
shading interp;
axis equal tight;
axis([0 2*pi 0 2*pi]);
ax = gca;ax.FontSize = 14;
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16);
title(['RCSX--iter2 2D, Time = ' num2str(Tf)], ...
    'Interpreter', 'latex', 'FontSize', 18);
colorbar;colormap parula;
% -------------------------------------------------------------
% Free Energy
% -------------------------------------------------------------
figure(2); clf;
set(gcf, 'Color', 'w');
idx = E_history > 0;
time_plot = time_axis(idx);
E_plot = E_history(idx);
loglog(time_plot, E_plot, 'LineWidth', 2.5);
grid on;grid minor;
ax = gca;ax.FontSize = 14;
ax.TickLabelInterpreter = 'latex';
xlabel('Time ($t$)', ...
    'Interpreter', 'latex', 'FontSize', 16);
ylabel('Free Energy $\mathcal{E}(u)$', ...
    'Interpreter', 'latex', 'FontSize', 16);
title(['2D Energy Dissipation | $N = ' num2str(N) '^2$'], ...
    'Interpreter', 'latex', 'FontSize', 18);
axis tight;
% -------------------------------------------------------------
% Final computing-time summary
% -------------------------------------------------------------
fprintf('\n---------------------------------------------\n');
fprintf('RCSX--iter2 2D Simulation Summary\n');
fprintf('---------------------------------------------\n');
fprintf('Grid:              %d^2\n', N);
fprintf('Time interval:     [0, %.2f]\n', Tf);
fprintf('Time step:         %.4f\n', dt);
fprintf('epsilon:           %.4f\n', eps);
fprintf('Computing time:    %.4f seconds\n', computing_time);
fprintf('Computing time:    %.4f minutes\n', minutes);
fprintf('Computing time:    %.4f hours\n', hours);
fprintf('Mass error:        %.16e\n', RealMassError);
fprintf('---------------------------------------------\n');