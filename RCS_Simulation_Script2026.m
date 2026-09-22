% Simulation Script to call the RCS 3D Solver
% For N=512, code will require a GPU with at least 24 Gbs of VRAM. 
clc
tic;
a=0;b=2*pi;
N = 256/2; Tf = 10; dt = 0.01; eps = 0.1;
[U, t_ax, E] = CH3D_RCSX_Solver(N, Tf, dt, eps);

% --- Domain setup ---
x = linspace(0, 2*pi, N);[X, Y, Z] = meshgrid(x, x, x);

rng(1527, 'twister');U0 = 0.05*rand(N,N,N) - 0.025;
% Direct initialization is used for the 3D implementation.
% The 50-step bootstrap is retained only in the 2D verification code.
InitialMass = sum(U0(:));

figure(1); %clf;
isosurface(X,Y,Z,U,-.3)
isosurface(X,Y,Z,U,-.15);
isosurface(X,Y,Z,U,-.05);
isosurface(X,Y,Z,U,.05);
isosurface(X,Y,Z,U,.15);
isosurface(X,Y,Z,U,.3)

ax = gca; ax.FontSize = 14;
camlight; lighting phong;
axis([0 2*pi 0 2*pi 0 2*pi]);
grid on; view(3);
title(['RCSX--iter2 3D, Time = ' num2str(Tf)]); 

% --- Gráfica de Energía 3D ---
figure(2); clf;
set(gcf, 'Color', 'w');

% Variables for energy computations: E y t_ax
idx = E > 0; time_plot = t_ax(idx);
E_plot = E(idx);

loglog(time_plot, E_plot, 'b-', 'LineWidth', 2.5);
grid on; grid minor;

ax = gca; ax.FontSize = 14;
ax.TickLabelInterpreter = 'latex';

xlabel('Time ($t$)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Free Energy $\mathcal{E}(u)$', 'Interpreter', 'latex', 'FontSize', 16);
title(['3D Energy Dissipation | $N = ' num2str(N) '^3$'], 'Interpreter', 'latex', 'FontSize', 18);
axis tight;
    % --- Finalizing Simulation ---
    U = gather(U); % Sending back to CPU
    FinalMass = sum(U(:));
    fprintf('\nSimulation Completed.\n');
    fprintf('Final Mass: %.16e\n', FinalMass);   
    % ---Priting/Verification of computations ---
h = (b-a)/N;Vol = h^3;RealMassError = abs(InitialMass - FinalMass) * Vol;
fprintf('Mass Error Computed: %.16e\n', RealMassError);    

%save('CH3D_Final_Result.mat', 'U', 'E_history', 'time_axis', 'InitialMass', 'FinalMass', 'MassError', 'N', 'dt', 'epsilon');
a=toc;minutes=a/60;hours=a/60^2;
minutes_hours=[minutes hours] 
