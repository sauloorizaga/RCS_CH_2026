% RCS 2D Solver for the Cahn--Hilliard Equation
% Self-contained GPU implementation
% RCSX-iter2 with 50-step bootstrap
clc
clear
close all
tic;
% --- Physical and numerical parameters ---
a=0;b=2*pi;
N=256;
Tf=100;
dt=0.01;
eps=0.05;
eps2=eps^2;
a_split=2;
h=(b-a)/N;
% --- GPU Mesh and Spectral Operators ---
x=gpuArray(linspace(a,b-h,N));
k=gpuArray([[0:N/2] [-N/2+1:-1]]);
[Kx,Ky]=meshgrid(k,k);
k2g=Kx.^2+Ky.^2;k4g=k2g.^2;
% --- Linear operators ---
lhs_cs=1+dt*(eps2*k4g+a_split*k2g);
lhs_f=1+(dt/2)*(eps2*k4g+a_split*k2g);
% --- Initial condition and bootstrap ---
rng(100,'twister');
Uinit = 0.05*rand(N, N) - 0.025; % Clean random start
Uinit=gpuArray(Uinit);
InitialMass=gather(sum(Uinit(:)));
m_boot=50;
dt_tiny=dt/m_boot;
lhs_tiny=1+dt_tiny*(eps2*k4g+a_split*k2g);
U_boot=Uinit;
for m=1:m_boot
    fU_boot=U_boot.^3-(1+a_split)*U_boot;
    U_boot=real(ifft2((fft2(U_boot)+dt_tiny*(-k2g.*fft2(fU_boot)))./lhs_tiny));
end
U=U_boot;
U_old=Uinit;
hat_U_n=fft2(U);
% --- Pre-allocation ---
steps=ceil(Tf/dt);
E_history=zeros(steps,1);
time_axis=zeros(steps,1);
fprintf('Starting 2D RCSX-iter2 (N=%d, dt=%.3f)...\n',N,dt);
t=0;it=1;
while t<Tf-dt*0.01
    U_extrap=2*U-U_old;
    % --- Richardson Solver ---
    Uc=U_extrap;
    for i=1:2
        fUc=Uc.^3-(1+a_split)*Uc;
        Uc=real(ifft2((hat_U_n+dt*(-k2g.*fft2(fUc)))./lhs_cs));
    end
    Uh=U_extrap;
    for i=1:2
        fUh=Uh.^3-(1+a_split)*Uh;
        Uh=real(ifft2((hat_U_n+(dt/2)*(-k2g.*fft2(fUh)))./lhs_f));
    end
    Uf=Uh;hat_Uh=fft2(Uh);
    for i=1:2
        fUf=Uf.^3-(1+a_split)*Uf;
        Uf=real(ifft2((hat_Uh+(dt/2)*(-k2g.*fft2(fUf)))./lhs_f));
    end
    U_next=2*Uf-Uc;
    % --- Energy ---
    hatU_next=fft2(U_next);
    E_grad=(eps2/2)*h^2*sum(k2g(:).*abs(hatU_next(:)).^2)/N^2;
    E_pot=0.25*h^2*sum((U_next(:).^2-1).^2);
    E_history(it)=gather(E_grad+E_pot);
    % --- Update ---
    U_old=U;
    U=U_next;
    hat_U_n=hatU_next;
    t=t+dt;
    time_axis(it)=t;it=it+1;
end
U=gather(U);
FinalMass=sum(U(:));
RealMassError=abs(InitialMass-FinalMass)*h^2;
fprintf('\nSimulation Completed.\n');
fprintf('Initial Mass: %.16e\n',InitialMass);
fprintf('Final Mass: %.16e\n',FinalMass);
fprintf('Mass Error Computed: %.16e\n',RealMassError);
time=toc;
fprintf('Runtime: %.2f min (%.2f hr)\n',time/60,time/3600);
% --- Plots ---
idx=E_history>0;
figure(1);
pcolor(gather(U));shading interp;axis equal tight;
ax=gca;ax.FontSize=14;
title(['RCSX--iter2 2D, Time = ' num2str(Tf)]);
figure(2);
loglog(time_axis(idx),E_history(idx),'b-','LineWidth',2.5);
grid on;grid minor;
ax=gca;ax.FontSize=14;ax.TickLabelInterpreter='latex';
xlabel('Time ($t$)','Interpreter','latex','FontSize',16);
ylabel('Free Energy $\mathcal{E}(u)$','Interpreter','latex','FontSize',16);
title(['2D Energy Dissipation | $N = ' num2str(N) '^2$'],'Interpreter','latex','FontSize',18);
axis tight;
