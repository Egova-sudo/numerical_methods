clc
clear all
close all

% Time dependent 2D Heat transfer problem 
% with homogeneous Drichlet boundary conditions
% Partial differential equation: T_t = T_xx + T_yy
% Use finite difference approximation of the second derivatives

spatial_step = [3,7,15,31];
time_step = [1/64,1/128,1/256,1/512,1/1024,1/2048,1/4096];
explicit_solution = cell(length(spatial_step),length(time_step));
implicit_solution = cell(length(spatial_step),length(time_step));

for m = 1:length(time_step)
    dt = time_step(m);
    time_domain = 0:dt:0.5;
    for n = 1:length(spatial_step)
        Nx = spatial_step(n); Ny = spatial_step(n);
        temperature_array_explicit = zeros(Nx+2,Ny+2,length(time_domain));
        temperature_array_implicit = zeros(Nx+2,Ny+2,length(time_domain));
        
        % Explicit Euler computation
        temperature_array_explicit(2:Nx+1,2:Ny+1,1) = 1; % initial condition
        for current_time_index = 1:length(time_domain)
            temperature_array_explicit(:,:,current_time_index+1) = ...
                explicitEuler(Nx,Ny,dt,temperature_array_explicit(:,:,current_time_index));
        end    
        
        % Implicit Euler computation
        % Skip computing by implicit euler step for time steps other than 1/64
        if time_step(m) == 1/64       
            temperature_array_implicit(2:Nx+1,2:Ny+1,1) = 1; % initial condition
            for current_time_index = 1:length(time_domain)
                temperature_array_implicit(:,:,current_time_index+1) = ...
                    implicitEuler(Nx,Ny,dt,temperature_array_implicit(:,:,current_time_index));
            end
        end
        
        % Saving results
        explicit_solution{n,m} = temperature_array_explicit;
        implicit_solution{n,m} = temperature_array_implicit;
    end
    
end


% Coloured surface map
time_checkpoints = [1/8,2/8,3/8,4/8];
methods = {explicit_solution,implicit_solution}; method_names = {'Explicit Euler','Implicit Euler'};
for method_index = 1:length(methods)
    method_name = method_names{method_index};
    method_solution = methods{method_index};
for plot_index = 1:length(time_checkpoints)
    % Creating a main figure window 
    fig = figure('Name',[method_name,' at t=',num2str(time_checkpoints(plot_index))]);
    if strcmp(method_name,'Implicit Euler')
        x0=10;y0=0;width=400;height=750;
        set(gcf,'position',[x0,y0,width,height])
    else
        x0=10;y0=0;width=1100;height=750;
        set(gcf,'position',[x0,y0,width,height])
    end
    position=1;
    for m = 1:length(spatial_step)
            Nx = spatial_step(m); Ny = spatial_step(m);
            hx = 1/(Nx+1); hy = 1/(Ny+1);
            xx = 0:hx:1; yy = 0:hy:1;
            for n = 1:length(time_step)
                % Skipping time steps other than 1/64 for Implicit Euler
                if strcmp(method_name,'Implicit Euler') && time_step(n) ~= 1/64
                    continue
                end
                dt = time_step(n); time_domain = 0:dt:0.5;
                % Creating a 4 by 7 window grip, each node will have a plot
                if strcmp(method_name,'Explicit Euler')
                    subplot(length(spatial_step),length(time_step),position)
                elseif strcmp(method_name,'Implicit Euler')
                    subplot(length(spatial_step),1,position)
                end
                position = position+1;
                time_checkpoint_index = find(time_domain == time_checkpoints(plot_index));
                surf(xx,yy,method_solution{m,n}(:,:,time_checkpoint_index))
                % Creating a seperate hidden only one plot and saving it
                hidden_fig = figure('visible','off');
                surf(xx,yy,method_solution{m,n}(:,:,time_checkpoint_index));
                saveas(hidden_fig,[method_name,' Nx_',int2str(Nx),' dt_',num2str(dt),...
                    ' t_',num2str(time_checkpoints(plot_index)),'.png']);
                close(hidden_fig)
            end
    end
end
end
