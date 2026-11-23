% =========================================================
% (3D): setup time, time per step and peak GPU memory
% BCP, RCSX-iter2 solver exactly as in BCP3D_RCSX_Solver_512
% Parameters of the 512^3 run reported in the paper (Fig. BCP_512)
% =========================================================
clc; clear;
N = 512; dt = 0.005; epsilon = 0.1; alpha = 38; ubar_offset = 0.34;
nsteps = 200;                         % short run: enough for time/step and peak memory

g = gpuDevice; reset(g); g = gpuDevice;   % clean device before measuring
memUsed = @() (g.TotalMemory - g.AvailableMemory)/1e9;   % GB used on the device
mem0 = memUsed();
fprintf('GPU: %s  (%.1f GB total)\n', g.Name, g.TotalMemory/1e9);

% ---------------- SETUP (timed) ----------------
wait(g); tSetup = tic;
eps2 = epsilon^2; a_split = 2; b_dom = 2*pi; h = b_dom/N;
k = gpuArray([[0:N/2] [-N/2+1:-1]]);
[Kx, Ky, Kz] = meshgrid(k, k, k);
k2g = Kx.^2 + Ky.^2 + Kz.^2;
inv_lhs_cs = 1 ./ (1 + dt*(eps2*(k2g.^2) + a_split*k2g + alpha));
inv_lhs_f  = 1 ./ (1 + (dt/2)*(eps2*(k2g.^2) + a_split*k2g + alpha));
clear Kx Ky Kz k;
rng(1527, 'twister');
U = gpuArray(0.01 * rand(N, N, N) + ubar_offset);
InitialSum = gather(sum(U(:)));
ubar = InitialSum / (N^3);
InitialSum = InitialSum*h^3;
hat_alpha_ubar_dt = zeros(N, N, N, 'gpuArray');
hat_alpha_ubar_dt(1,1,1) = alpha * ubar * (N^3) * dt;
hat_U_n = fftn(U); U_old = U;
wait(g); setupTime = toc(tSetup);
memSetup = memUsed();

% ---------------- TIME STEPPING (timed) ----------------
peakMem = memSetup; stepTimes = zeros(nsteps,1);
for it = 1:nsteps
    wait(g); ts = tic;
    U_extrap = 2*U - U_old;
    Uc = U_extrap;
    for i = 1:2
        Uc = real(ifftn((hat_U_n + dt*(-k2g.*fftn(2*(Uc.^3 - Uc) - a_split*Uc)) + hat_alpha_ubar_dt) .* inv_lhs_cs));
    end
    Uh = U_extrap;
    for i = 1:2
        Uh = real(ifftn((hat_U_n + (dt/2)*(-k2g.*fftn(2*(Uh.^3 - Uh) - a_split*Uh)) + 0.5*hat_alpha_ubar_dt) .* inv_lhs_f));
    end
    Uf = Uh; hat_Uh = fftn(Uh);
    for i = 1:2
        Uf = real(ifftn((hat_Uh + (dt/2)*(-k2g.*fftn(2*(Uf.^3 - Uf) - a_split*Uf)) + 0.5*hat_alpha_ubar_dt) .* inv_lhs_f));
    end
    U_old = U;
    U = 2*Uf - Uc;
    hat_U_n = fftn(U);
    wait(g); stepTimes(it) = toc(ts);
    peakMem = max(peakMem, memUsed());
end

FinalSum = gather(sum(U(:)))*h^3;
tFirst = stepTimes(1);
tStep  = mean(stepTimes(11:end));     % steady-state time per step (excludes warm-up)
estTotal = setupTime + tFirst + (20/dt - 1)*tStep;

fprintf('\n================ RESULTS (N = %d^3) ================\n', N);
fprintf('Setup time (grid, operators, IC on GPU): %.2f s\n', setupTime);
fprintf('First step (includes GPU/FFT warm-up):   %.2f s\n', tFirst);
fprintf('Time per step (steady state):           %.3f s\n', tStep);
fprintf('Estimated total to T_f = 20 (%d steps):  %.1f min\n', round(20/dt), estTotal/60);
fprintf('GPU memory before setup:                %.2f GB\n', mem0);
fprintf('GPU memory after setup:                 %.2f GB\n', memSetup);
fprintf('Peak GPU memory during stepping:        %.2f GB\n', peakMem);
fprintf('Relative mass error after %d steps:     %.2e\n', nsteps, abs(FinalSum-InitialSum)/abs(InitialSum));
fprintf('=====================================================\n');
