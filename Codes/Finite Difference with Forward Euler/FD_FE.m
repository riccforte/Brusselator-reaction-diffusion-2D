%% Project CCRE-ATP_APC
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

% Visualize grid
figure(1)
plot(X,Y,'.')
title('Meshgrid')

%% Problem Data
%Input Parameters
A = 3.4;
B = 1;           %[-]
alpha = 0.002;     %[-]    

dt = 0.01;       %[s]     %time step
tau = 5;         %[s]     %total time of simulation
nsteps = tau/dt;  %[-]     %total time steps

%% Initial conditions
u = 0.5 + Y;
v = 1 + 5 .* X;

figure(2)
surf(X,Y,u)
xlabel('X')
ylabel('Y')
zlabel('u')
title('U initial conditions')

figure(3)
surf(X,Y,v)
surf(X,Y,u)
xlabel('X')
ylabel('Y')
zlabel('v')
title('V initial conditions')


%% Time iterative loop
t = 0;
u_storage = zeros(N, M, nsteps);
v_storage = zeros(N, M, nsteps);

figure
for index = 1 : nsteps
    uo = u;     
    vo = v;

%--------------------------------------------------------------------------|
    for i = 2 : N-1
        for j = 2:M-1  
%---------------------------- Solving for u--------------------------------|
         dudx2 = (u(i+1,j) - 2*uo(i,j) + u(i-1,j)) / ((x(i+1) - x(i)))^2;
         dudy2 = (u(i,j+1) - 2*uo(i,j) + u(i,j-1)) / ((y(j+1) - y(j)))^2;
         deriv = dudx2 + dudy2;      
         u(i,j) = uo(i,j) + ( B + uo(i,j)^2*vo(i,j) - (A+1)*uo(i,j) + alpha*deriv )*dt;
%---------------------------- Solving for v--------------------------------|      
         dvdx2 = (v(i+1,j) - 2*vo(i,j) + v(i-1,j)) / (x(i+1) - x(i))^2;
         dvdy2 = (v(i,j+1) - 2*vo(i,j) + v(i,j-1)) / (y(j+1) - y(j))^2;
         deriv_v = dvdx2 + dvdy2; 
         v(i,j) = vo(i,j) + ( A*u(i,j) - u(i,j)^2*v(i,j) + alpha*deriv_v ) *dt;
        end
    end

  % Apply boundary conditions for u
    u(1,:) = u(2,:);  
    u(N,:) = u(N-1,:); 
    u(:,1) = u(:,2); 
    u(:,M) = u(:,M-1);   
 % Apply boundary conditions for v
   v(1,:) = v(2,:);  
   v(N,:) = v(N-1,:); 
   v(:,1) = v(:,2); 
   v(:,M) = v(:,M-1);
   
  %Store values
  u_storage(:, :, index) = u;
  v_storage(:, :, index) = v;
  
  %Update time
  t = t + dt;
  %Live plot
  subplot(2,1,1)
  surf(X,Y,u);
  title(['u at Time: ' num2str(t) 's']);
  subplot(2,1,2)
  surf(X,Y,v);
  title(['v at Time: ' num2str(t) 's']);
  drawnow;
  disp(num2str(t));
end
 %------------------------------------------------------------------------%
  
    
   




