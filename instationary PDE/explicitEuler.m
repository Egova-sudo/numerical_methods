function result = explicitEuler(Nx,Ny,dt,T)
hx = 1/(Nx+1);
hy = 1/(Ny+1);
result = zeros(Nx+2,Ny+2);

new_T = T(2:Nx+1,2:Ny+1) + dt * (T(1:Nx,2:Ny+1) -2*T(2:Nx+1,2:Ny+1) ...
    +T(3:Nx+2,2:Ny+1))/hx^2 + dt*(T(2:Nx+1,1:Ny) -2*T(2:Nx+1,2:Ny+1) ...
    +T(2:Nx+1,3:Ny+2))/hy^2;
result(2:Nx+1,2:Ny+1) = new_T;
end