%% Project CCRE-ATP-APC
%Clearing Matlab Environment
clear
clc
close all

%% Grid Generation

%2D Domain
N = 20;    %Points along x direction
M = 20;    %Points along y direction

%Length Domain
Lx = 1;
Ly = 1;

%Grid Transformation
xi = @(i) 0.5*(1-cos((i-1)*pi/(N-1)))*Lx;
yi = @(j) 0.5*(1-cos((j-1)*pi/(M-1)))*Ly;

%Preallocating Final Grid Variables
x = zeros(N,1);
y = zeros(M,1);

%Calculating grid points
for i = 1: N
   x(i) = xi(i);
end
for j = 1:M
    y(j) = yi(j);
end

%Creating rectangular Grid
[X ,Y] =meshgrid(x,y);

% % Visualize grid
% figure(1)
% plot(X,Y,'.')
% title('Meshgrid')

%% Problem Data
%Input Parameters
A = 3.4;
B = 1;             %[-]
alpha = 0.002;       %[-]    

dt = 0.01;        %[s]     %time step
tau = 10;           %[s]     %total time of simulation
nsteps = tau/dt;   %[-]     %total time steps

%% Initial conditions

%Preallocating Variables
u = zeros(N,M);
v = zeros(N,M);
u_storage = zeros(N, M, nsteps);
v_storage = zeros(N, M, nsteps);

for i = 1:N
    for j = 1:M
        u(i,j) = 0.5 + y(j);
        v(i,j) = 1 + 5*x(i);
    end
end
figure
subplot(2,1,1)
surf(X,Y,u)
subplot(2,1,2)
surf(X,Y,v)

%% Calculating the weighting coefficient for the 1st and 2nd derivative

% X direction
[a1] = weighted_coefficients1(N,M,x);
[a2] = weighted_coefficient2(N,M,a1,x);

% Y direction
[b1] = weighted_coefficients1(N,M,y);
[b2] = weighted_coefficient2(N,M,b1,y);

%% Functions
t = 0;
h = figure;
  
% Main simulation loop using RK4
u_old = u;  % Initialize u_old with the initial condition of u
v_old = v;  % Initialize v_old with the initial condition of v

for time = 1 : nsteps 
    % Apply boundary conditions at the start of the time step if necessary
    [u_old, v_old] = updateMatrices(u_old, v_old, a1, b1);
    
    % RK4 stage 1
    [k1_u, k1_v] = computeRHS(u_old, v_old, a1, a2, b1, b2, alpha, A, B, N, M);
    
    % RK4 stage 2
    [k2_u, k2_v] = computeRHS(u_old + 0.5*dt*k1_u, v_old + 0.5*dt*k1_v, a1, a2, b1, b2, alpha, A, B, N, M);
    
    % RK4 stage 3
    [k3_u, k3_v] = computeRHS(u_old + 0.5*dt*k2_u, v_old + 0.5*dt*k2_v, a1, a2, b1, b2, alpha, A, B, N, M);
    
    % RK4 stage 4
    [k4_u, k4_v] = computeRHS(u_old + dt*k3_u, v_old + dt*k3_v, a1, a2, b1, b2, alpha, A, B, N, M);
    
    % Update the u and v with the new values
    u = u_old + (dt/6)*(k1_u + 2*k2_u + 2*k3_u + k4_u);
    v = v_old + (dt/6)*(k1_v + 2*k2_v + 2*k3_v + k4_v);

    % Apply boundary conditions again if they are dynamic and change with each time step
    [u, v] = updateMatrices(u, v, a1, b1);
    
    % Prepare for the next time step
    u_old = u;
    v_old = v;
     
         %Live plot
         %u
         set(h, 'Position', [280 150 1000 600]);
         subplot(2,2,1)
         surf(X,Y,u)
         title(['U surface at ' num2str(t)])
         zlabel('V')
         xlabel('X')
         ylabel('Y')
         %v
         subplot(2,2,2)
         surf(X,Y,v)
         title(['V surface at ' num2str(t)])
         zlabel('V')
         xlabel('X')
         ylabel('Y')
         drawnow
         subplot(2,2,3)
         contourf(X,Y,u)
         subplot(2,2,4)
         contourf(X,Y,v)
                
%     %Storing every value of u
%     u_storage(:, :, time) = u;
%     v_storage(:, :, time) = v;
    %Time counter
    t = t + dt; 
end

%% Functions
function [RHS_u, RHS_v] = computeRHS(u, v, a1, a2, b1, b2, alpha, A, B, N, M)

    RHS_u = zeros(N,M);
    RHS_v = zeros(N,M);
    Derivatives_u = zeros(N,M);
    Derivatives_v = zeros(N,M);
    
    % Computation for RHS of u
    for i = 2 : N-1
        for j = 2:M-1
            Derivatives_u(i,j) = sum(a2(i,:) .* u(:,j)') + sum(b2(j,:) .* u(i,:));
            RHS_u(i,j) = B + u(i,j)^2*v(i,j) - (A+1)*u(i,j) + alpha*Derivatives_u(i,j);
           
            Derivatives_v(i,j) = sum(a2(i,:) .* v(:,j)') + sum(b2(j,:) .* v(i,:));
            RHS_v(i,j) = A*u(i,j) - u(i,j)^2*v(i,j) + alpha*Derivatives_v(i,j);
        end
    end

end
