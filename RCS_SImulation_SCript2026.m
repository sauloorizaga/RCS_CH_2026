% --- Run Test ---
tic;
a=0;b=2*pi;
N = 256/2; Tf = 10; dt = 0.01; eps = 0.1;
[U, t_ax, E] = CH3D_RCSX_Solver(N, Tf, dt, eps);

% --- Plotting Morfología ---
x = linspace(0, 2*pi, N); 
[X, Y, Z] = meshgrid(x, x, x);

rng(1527, 'twister');
U0 = 0.05*rand(N,N,N) - 0.025;
InitialMass = sum(U0(:));

figure(3); clf;
%set(gca, 'NextPlot', 'add'); % El reemplazo elegante del hold on

isosurface(X,Y,Z,U,-.3)
isosurface(X,Y,Z,U,-.15);
isosurface(X,Y,Z,U,-.05);
isosurface(X,Y,Z,U,.05);
isosurface(X,Y,Z,U,.15);
isosurface(X,Y,Z,U,.3)

ax = gca; 
ax.FontSize = 14;
camlight; lighting phong;
axis([0 2*pi 0 2*pi 0 2*pi]);
grid on; view(3);
title(['RCSX--iter2 3D, Time = ' num2str(Tf)]); 

% --- Gráfica de Energía 3D ---
figure(600); clf;
set(gcf, 'Color', 'w');

% USAMOS LAS VARIABLES QUE ESCUPE LA FUNCIÓN: E y t_ax
idx = E > 0; 
time_plot = t_ax(idx);
E_plot = E(idx);

loglog(time_plot, E_plot, 'b-', 'LineWidth', 2.5);
grid on; grid minor;

ax = gca; 
ax.FontSize = 14;
ax.TickLabelInterpreter = 'latex';

xlabel('Time ($t$)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Free Energy $\mathcal{E}(u)$', 'Interpreter', 'latex', 'FontSize', 16);
title(['3D Energy Dissipation | $N = ' num2str(N) '^3$'], 'Interpreter', 'latex', 'FontSize', 18);
axis tight;

%MASS

    % --- Finalizing Simulation ---
    U = gather(U); % La bajamos al CPU para que el .mat no pese tanto y sea compatible
    FinalMass = sum(U(:));
     fprintf('\nSimulacion Terminada.\n');
    fprintf('Masa Final: %.16e\n', FinalMass);
    
    % --- En tu script de finalización ---
h = (b-a)/N;
Vol = h^3;
RealMassError = abs(InitialMass - FinalMass) * Vol;
fprintf('Error de Masa Real: %.16e\n', RealMassError);

    
 
    
    % EL SALVAVIDAS: Guarda todo en un archivo con el timestamp o nombre fijo
    %save('CH3D_Final_Result.mat', 'U', 'E_history', 'time_axis', 'InitialMass', 'FinalMass', 'MassError', 'N', 'dt', 'epsilon');
   
%-------------

a=toc;
minutes=a/60;
hours=a/60^2;
minutes_hours=[minutes hours]
 
