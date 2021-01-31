function result = implicitEuler(Nx,Ny,dt,T)
hx = 1/(Nx+1); hy = 1/(Ny+1); N = Nx*Ny;
k1 = 1/(1 + 2*dt/hx^2 + 2*dt/hy^2);
k2 = dt/hx^2; k3 = dt/hy^2;
new_T = zeros(Nx+2,Ny+2);

iteration = 0; error = 1; k4_old = 0;
% Gauss-Seidel iteration
while error >= 1e-6
    iteration = iteration +1;
    % Implicit Euler part
    new_T(2:Nx+1,2:Ny+1) = k1 * (T(2:Nx+1,2:Ny+1) + ...
        k2 * (new_T(1:Nx,2:Ny+1) + new_T(3:Nx+2,2:Ny+1))+ ...
        k3 * (new_T(2:Nx+1,1:Ny) + new_T(2:Nx+1,3:Ny+2))) ;
    
    % Calculation residual norm as error
    k4 = T(2:Nx+1,2:Ny+1) +k2 * (new_T(1:Nx,2:Ny+1) -2*new_T(2:Nx+1,2:Ny+1) ...
    +new_T(3:Nx+2,2:Ny+1)) + k3*(new_T(2:Nx+1,1:Ny) -2*new_T(2:Nx+1,2:Ny+1) ...
    +new_T(2:Nx+1,3:Ny+2));
    S = sum((k4-k4_old).^2,'all');
    error = sqrt(S/N);
    k4_old = k4;
end

result = new_T;
end