function [U, time_axis, E_history] = CH3D_RCSX_Solver(N, tfinal, dt, epsilon)
    % --- Physical Parameters ---
    eps2 = epsilon^2;
    a_split = 1; 
    M = 1; a = 0; b = 2*pi; h = (b-a)/N;
    
    % --- GPU Mesh and Spectral Operators ---
    x = gpuArray(linspace(a, b-h, N));
    % [X, Y, Z] = meshgrid(x, x, x); % Borrado: no se ocupa dentro del loop (ahorra VRAM)
    
    k = gpuArray([[0:N/2] [-N/2+1:-1]] ./ M);
    [Kx, Ky, Kz] = meshgrid(k, k, k);
    k2g = Kx.^2 + Ky.^2 + Kz.^2;
    k4g = k2g.^2;
    
    lhs_cs = 1 + dt*(eps2*k4g + a_split*k2g);
    lhs_f  = 1 + (dt/2)*(eps2*k4g + a_split*k2g);
    
    rng(1527, 'twister');
    U = gpuArray(0.05*rand(N, N, N) - 0.025);
    hat_U_n = fftn(U);
    U_old = U; 
    
    
    steps = ceil(tfinal/dt);
    E_history = zeros(steps, 1);
    time_axis = zeros(steps, 1);
    
    fprintf('Starting 3D RCSX-iter2 (N=%d, dt=%.3f)...\n', N, dt);
    tic;
        t = 0; it = 1;   
        % Grab initial mass (Corrected line)
InitialMass = gather(sum(U(:))); 
fprintf('Masa Inicial: %.16e\n', InitialMass)
    
    while t < tfinal - dt*0.01
        U_extrap = 2*U - U_old; 
        
        % --- Richardson Solver ---
        Uc = U_extrap;
        for i = 1:2
            fUc = Uc.^3 - (1+a_split)*Uc;
            Uc = real(ifftn((hat_U_n + dt*(-k2g.*fftn(fUc))) ./ lhs_cs));
        end
        
        Uh = U_extrap;
        for i = 1:2
            fUh = Uh.^3 - (1+a_split)*Uh;
            Uh = real(ifftn((hat_U_n + (dt/2)*(-k2g.*fftn(fUh))) ./ lhs_f));
        end
        
        Uf = Uh; hat_Uh = fftn(Uh);
        for i = 1:2
            fUf = Uf.^3 - (1+a_split)*Uf;
            Uf = real(ifftn((hat_Uh + (dt/2)*(-k2g.*fftn(fUf))) ./ lhs_f));
        end
        
        U_next = 2*Uf - Uc;
        
        % --- ENERGY (Parseval's version: Fast & Precise) ---
        if mod(it, 1) == 0 || it == 1
            hatU_next = fftn(U_next); 
            E_grad = (eps2 / 2) * (h^3) * sum(k2g(:) .* abs(hatU_next(:)).^2) / (N^3);
            E_pot = 0.25 * (h^3) * sum((U_next(:).^2 - 1).^2);
            E_history(it) = gather(E_grad + E_pot);
        end
        
        % --- Update ---
        U_old = U;
        U = U_next;
        hat_U_n = fftn(U);
        t = t + dt;
        time_axis(it) = t;
        it = it + 1;
    end
    U = gather(U); 
    toc;
end
