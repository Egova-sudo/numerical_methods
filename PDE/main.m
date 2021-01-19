clc
clear all
close all

% Two-dimensional stationary heat equation:
% T_xx + T_yy = -2(pi)**2 sin(pi*x) sin(pi*y)
% Boundary condition:
% T(x,y) = 0 for edges

% Parameters
steps = [3,7,15,31,63,127];
T_solutions =cell(4); runtime = zeros(length(steps),3); storage = zeros(length(steps),3);
for k = 1:length(steps)
N_x = steps(k); N_y = steps(k);
N = N_x*N_y;

% Right hand side, b
bb = zeros(N_x+2,N_y+2);
func = @(x,y) -2*pi^2 * sin(pi*x)' .* sin(pi*y);
h_x = 1/(N_x+1); h_y= 1/(N_y+1);
xx = h_x:h_x:h_x*N_x; yy = h_y:h_y:h_y*N_y;
bb_dummy = func(xx,yy);
bb(2:N_x+1,2:N_y+1) = bb_dummy;
b_gauss = reshape(bb,[],1);
b = reshape(bb_dummy,[],1);

% Analytical solution:
f_analytical = @(x,y) sin(pi*x)' .* sin(pi*y);
T_analytical = f_analytical(xx,yy);

% Creating a matrix from linear system of equations
A = LinearEquationSystem(N_x,N_y);

% Full matrix direct solution
ti_full = tic;
T_full = A\b;
T_full = reshape(T_full,N_y,N_x)';
tf_full = toc(ti_full);

% Sparse matrix direct solution
ti_sparse = tic;
A_sparse = sparse(A);
T_sparse = A_sparse\b;
T_sparse = reshape(T_sparse,N_y,N_x)';
tf_sparse = toc(ti_sparse);

% Gauss-seidel Iterative solver
ti_gauss = tic;
T_gauss = GaussSeidel(b_gauss,N_x,N_y);
tf_gauss = toc(ti_gauss);

% Visualization of the results
if N_x~=127
xi = 0:h_x:1; yi= 0:h_y:1;
methods = {T_full, T_sparse, T_gauss};
methods_name = {'Full Matrix with direct sol','Sparse matrix with direct sol','Gauss-Seidel iterative sol'};

figure 
for i=1:length(methods)
    solution = zeros(N_x+2,N_y+2); % adding boundary condition
    solution(2:N_x+1,2:N_y+1) = methods{i};
   
    % Coloured surface map
    subplot(length(methods),2,2*i-1)
    surf(xi,yi,solution);
    title(strcat(methods_name{i},' for N_x=',int2str(N_x)));
    % Countor map
    subplot(length(methods),2,2*i)
    contour(xi,yi,solution);
    title(strcat(methods_name{i},' for N_x=',int2str(N_x)));
end
end
% End of loop
T_solutions{k,1}=T_full; T_solutions{k,2}=T_sparse; T_solutions{k,3}=T_gauss; T_solutions{k,4} =T_analytical;
runtime(k,:) = [tf_full,tf_sparse,tf_gauss];
storage_A = whos('A').bytes/8; % that gives the number of doubles in the matrix A
storage_A_sparse = whos('A_sparse').bytes/8; 
storage_T_gauss = whos('T_gauss').bytes/8;
storage_b = whos('b').bytes/8;
storage_b_gauss = whos('b_gauss').bytes/8;
storage(k,:) = [storage_A + storage_b, storage_A_sparse+storage_b, storage_T_gauss+storage_b_gauss];
end

% Runtime and storage comparison for N_x = 7,15,31,and 63
table([7,15,31,63;runtime(2:5,1)';storage(2:5,1)'],'RowNames',{'N_x,N_y';'runtime';'storage'},'VariableNames',{'direct solution with full matrix'})
table([7,15,31,63;runtime(2:5,2)';storage(2:5,2)'],'RowNames',{'N_x,N_y';'runtime';'storage'},'VariableNames',{'direct solution with sparse matrix'})
table([7,15,31,63;runtime(2:5,3)';storage(2:5,3)'],'RowNames',{'N_x,N_y';'runtime';'storage'},'VariableNames',{'iterative solution with Gauss-seidel'})

% Compute Gauss-Seidel errors and compare with the others
error = zeros(1,length(steps));
for i = 1:length(steps)
    N = steps(i)^2;
    total = sum((T_solutions{i,3} - T_solutions{i,4}).^2,'all');
    error(1,i) = sqrt(total/N);
end
error_reduction = [NaN,error(1:end-1)./error(2:end)];
% Tabularise the results
table([7,15,31,63,127;error(2:end);error_reduction(2:end)],'RowNames',{'N_x,N_y';'error';'error reduction'},'VariableNames',{'Gauss-Seidel'})
