clc
clear all
close all

global time_domain t_start
% Examining an ordinary differential equation defining the population 
% dynamic of a corresponding specie

f = @(p) (1 - p./10).*p; % example ODE
f_analytical = @(t) 10 ./ (1 + 9.*exp(-t)); % analytical solution of ODE

% Visualizion of how the function looks on time domain
t_start = 0; t_end = 5; 
time_domain = t_start:1/8:t_end;
y_exact = f_analytical(time_domain);
figure
plot(time_domain, y_exact);
ylabel('p(t)')
xlabel('time (t)')

% ODE is numerically solved and plotted on methods
% of Euler, Heun, 4th order Runge-kutta

% initial conditions
dt=[1,1/2,1/4,1/8]; y0 = 1;
method_name = ["Euler","Heun","RK"];
methods = {@euler,@heun,@rk}; % creating a cell array containing
                              % function handles of euler, heun, and rk
                              % methods
app_error = zeros(length(method_name),length(dt));

for i=1:length(method_name)
    numerical_method = methods{i};
    
    figure
    plot(time_domain,y_exact);
    for j=1:length(dt)
        hold on
        y = numerical_method(f, y0, dt(j), t_end);
        app_error(i,j) = approximation_error(y,y_exact,dt(j),t_end);
        plot(t_start:dt(j):t_end, y);
    end
    title(method_name(i)+' method')
    legend('\deltat_1='+string(dt(1)),'\deltat_2='+string(dt(2)),'\deltat_3='+string(dt(3)),'\deltat_4='+string(dt(4)),'Analytical','Location','southeast')
    saveas(gcf,method_name(i)+'.pdf')
    
    % Create a table for error comparison of step size
    table_error(i,[dt;app_error(i,:)]);
end
%% Functions

% Consider a general initial value problem
% y' = f(y); y(0) = y_0:
% Implement the following explicit numerical methods
% 1) explicit Euler method,
% 2) method of Heun,
% 3) Runge-Kutta method (fourth order)
% as a matlab function depending on the right hand side f(y)^1, the initial value y_0,
% the timestep size \delt and the end time t_end. The output of the function shall be a
% vector containing all computed approximate values for y.


function y = euler(y_prime, y0, dt, t_end)
y = zeros(t_end / dt +1,1);
t = dt; index = 1; y(1) = y0; % initial conditions
while t <= t_end
    y(index+1) = y(index) + dt*y_prime(y(index));
    index = index+1;t=t+dt;
end
end

function y = heun(y_prime, y0, dt, t_end)
y = zeros(t_end / dt +1,1);
t = dt; index = 1; y(1) = y0; % initial conditions
while t <= t_end
    k_1 = y_prime(y(index));
    y(index+1) = y(index) + dt*(1/2)*(k_1 + y_prime(y(index)+dt*k_1));
    index =index+1;t=t+dt;
end
end

function y = rk(y_prime, y0, dt, t_end)
y = zeros(t_end / dt +1,1);
t = dt; index = 1; y(1) = y0; % initial conditions
while t <= t_end
    k_1 = y_prime(y(index));
    k_2 = y_prime(y(index)+0.5*dt*k_1);
    k_3 = y_prime(y(index)+0.5*dt*k_2);
    k_4 = y_prime(y(index)+k_3*dt);
    y(index+1) = y(index) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt;  % main equation
    index =index+1;t=t+dt;
end
end

% Error computations
% Approximation errors
function result = approximation_error(p,p_exact,dt,t_end)
global time_domain t_start
% check the shape of p and p_exact. If it does not have a vector shape (1,~)
% Then take its transpose
if size(p,2) == 1
    p = p';
elseif size(p_exact,2) == 1
    p_exact = p_exact';
end
% size of p depends on the dt selection and size of p_exact
% comes from fplot command. These sizes may not match everytime.
% Thus, it should be handled. Think of better way?
% Looking for y values with the same t_k.
time_p = t_start:dt:t_end;
matching_index = ismember(time_domain,time_p);
difference = (p-p_exact(matching_index)).^2;
result = sqrt((dt/t_end)*sum(difference));

end

function T = table_error(i,data)
% Options for nice visualization
titles = ["Euler's explicit method O(h)","Heun's method O(h^2)","Runge-Kutta method O(h^4)"];
row = ["\deltat";"Error";"Error reduction"];
row = cellstr(row);

% Evaluating the error reduction by step size
[s1,s2] = size(data);
reduction = data(2,1:s2-1)./data(2,2:s2);
data = [data;nan,reduction];
T = table(data,'VariableNames',titles(i),'RowNames',row)
end