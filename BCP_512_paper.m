%% BCP_512_paper.m
% Richardson--Convex Splitting (RCS) solver for the 3D Block Copolymer
% (Ohta--Kawasaki) model on the GPU. Self-contained script.
%
% Default parameters reproduce Figure 14 of the paper (512^3 spherical
% microphase separation). Users can also explore the BCP model presented
% in the paper by changing the parameters below, e.g. the gyroid of
% Figures 12-13: N = 256, Tf = 100, alpha = 20, ubar_offset = 0.20.
%
% Requires MATLAB with Parallel Computing Toolbox and an NVIDIA GPU.
% At N = 512 (double precision) the run needs about 22 GB of GPU memory;
% use N = 128 or 256 on smaller GPUs.
%
% S. Orizaga, "A Second-Order Richardson--Convex Splitting Method for the
% Cahn--Hilliard Equation: Stability Analysis and GPU-Accelerated 3D
% Computations".

clear; clc;

%% Parameters
N           = 512/2;     % grid points per dimension
Tf          = 20;      % final time
dt          = 0.01;    % time step
epsilon     = 0.1;     % interface parameter
alpha       = 38;      % nonlocal (long-range) parameter
ubar_offset = 0.34;    % mean composition
a_split     = 1;       % convex splitting parameter

doEnergy    = true;    % compute the energy during the run
nE          = 10;      % compute the energy every nE steps

%% Setup
g = gpuDevice;
tic;
eps2 = epsilon^2; L = 2*pi; h = L/N;

% Wavenumbers
k = gpuArray([0:N/2, -N/2+1:-1]);
[Kx, Ky, Kz] = meshgrid(k, k, k);
k2g = Kx.^2 + Ky.^2 + Kz.^2;
clear Kx Ky Kz;                       % free memory

% Inverse of the implicit operators (coarse step and half-steps)
inv_lhs_cs = 1 ./ (1 + dt    *(eps2*k2g.^2 + a_split*k2g + alpha));
inv_lhs_f  = 1 ./ (1 + (dt/2)*(eps2*k2g.^2 + a_split*k2g + alpha));

% Initial condition: small random perturbation of the mean composition
rng(1527, 'twister');
U = gpuArray(0.01*rand(N, N, N) + ubar_offset);

InitialSum = gather(sum(U(:)));
ubar = InitialSum/N^3;
InitialMass = InitialSum*h^3;

% Zero-mode forcing alpha*ubar (mass term)
hat_alpha_ubar_dt = zeros(N, N, N, 'gpuArray');
hat_alpha_ubar_dt(1,1,1) = alpha*ubar*N^3*dt;

hat_U_n = fftn(U); U_old = U;
wait(g); tSetup = toc;
fprintf('Setup time: %.2f s\n', tSetup);
fprintf('Initial mass: %.16e\n', InitialMass);

%% Energy:  E = int( eps^2/2 |grad u|^2 + F(u) + alpha/2 |grad psi|^2 ),
%  F(u) = (u^2-1)^2/2,  -Lap psi = u - ubar  (scaling of the implemented equation)
nSteps = round(Tf/dt);
nRec   = floor(nSteps/nE) + 1;
t_ax = zeros(nRec, 1); E = zeros(nRec, 1); iRec = 1;
if doEnergy
    E(1) = bcpEnergy(U, hat_U_n, k2g, eps2, alpha, h, N);
end

%% Time stepping
tic;
for n = 1:nSteps
    U_extrap = 2*U - U_old;           % predictor

    % 1. Coarse step (dt), two fixed-point sweeps
    Uc = U_extrap;
    for i = 1:2
        Uc = real(ifftn((hat_U_n + dt*(-k2g.*fftn(2*(Uc.^3 - Uc) - a_split*Uc)) ...
             + hat_alpha_ubar_dt) .* inv_lhs_cs));
    end

    % 2. First half-step (dt/2)
    Uh = U_extrap;
    for i = 1:2
        Uh = real(ifftn((hat_U_n + (dt/2)*(-k2g.*fftn(2*(Uh.^3 - Uh) - a_split*Uh)) ...
             + 0.5*hat_alpha_ubar_dt) .* inv_lhs_f));
    end

    % 3. Second half-step (dt/2)
    Uf = Uh; hat_Uh = fftn(Uh);
    for i = 1:2
        Uf = real(ifftn((hat_Uh + (dt/2)*(-k2g.*fftn(2*(Uf.^3 - Uf) - a_split*Uf)) ...
             + 0.5*hat_alpha_ubar_dt) .* inv_lhs_f));
    end

    % Richardson combination
    U_old = U;
    U = 2*Uf - Uc;
    hat_U_n = fftn(U);

    if doEnergy && mod(n, nE) == 0
        clear Uc Uh Uf hat_Uh U_extrap;   % free memory before the energy
        iRec = iRec + 1;
        t_ax(iRec) = n*dt;
        E(iRec) = bcpEnergy(U, hat_U_n, k2g, eps2, alpha, h, N);
    end

    if mod(n, 100) == 0
        wait(g);
        fprintf('Step %5d / %d   t = %6.2f   elapsed %.1f min\n', n, nSteps, n*dt, toc/60);
    end
end
wait(g); tSolve = toc;

U = gather(U);
FinalMass = sum(U(:))*h^3;
fprintf('\n==========================================\n');
fprintf('  RCS BCP run complete (N = %d)\n', N);
fprintf('  Solver time:  %.2f min (%.3f s/step)\n', tSolve/60, tSolve/nSteps);
fprintf('  Final mass:   %.16e\n', FinalMass);
fprintf('  Mass error:   %.3e\n', abs(InitialMass - FinalMass));
fprintf('==========================================\n');

%% Visualization (isosurface rendering at 512^3 can take several minutes)
x = linspace(0, L, N);
figure(1); clf;
isosurface(x, x, x, U, -0.15);
isosurface(x, x, x, U, -0.05);
isosurface(x, x, x, U,  0.05);
isosurface(x, x, x, U,  0.15);
camlight; lighting phong;
axis([0 L 0 L 0 L]); set(gca, 'FontSize', 14);
title(sprintf('BCP, RCS, N = %d^3, T = %g', N, Tf));

if doEnergy
    t_ax = t_ax(1:iRec); E = E(1:iRec);
    figure(2); clf;
    semilogx(t_ax(2:end), E(2:end), 'b-', 'LineWidth', 2.5); grid on;
    xlabel('Time ($t$)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Energy $\mathcal{E}(u)$', 'Interpreter', 'latex', 'FontSize', 14);
    title('Energy Dissipation (Ohta--Kawasaki)', 'Interpreter', 'latex', 'FontSize', 16);
end

%% Local function
function En = bcpEnergy(U, hat_U, k2g, eps2, alpha, h, N)
    P = abs(hat_U).^2; P(1,1,1) = 0;              % remove the mean mode
    Egrad = 0.5*eps2*gather(sum(P.*k2g, 'all'));  % eps^2/2 |grad u|^2
    k2g(1,1,1) = 1;                               % local copy, avoid 0/0
    Enl   = 0.5*alpha*gather(sum(P./k2g, 'all')); % alpha/2 |grad psi|^2
    Ewell = 0.5*gather(sum((U.^2 - 1).^2, 'all'));
    En = h^3*(Ewell + (Egrad + Enl)/N^3);
end
