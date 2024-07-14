%% Project CCRE-ATP-APC
%Clearing Matlab Environment
clear
clc
close all

%% Grid Generation

%2D Domain
N = 25;    %Points along x direction
M = 25;    %Points along y direction

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
tau = 6;           %[s]     %total time of simulation
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

% figure(2)
% surf(X,Y,u)
% xlabel('X')
% ylabel('Y')
% zlabel('u')
% title('U initial conditions')
% 
% figure(3)
% surf(X,Y,v)
% xlabel('X')
% ylabel('Y')
% zlabel('v')
% title('V initial conditions')

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
Derivatives_u = zeros(N,M);
Derivatives_v = zeros(N,M); 

  for time = 1 : nsteps 
         
    % Computation for RHS of u
    for i = 2 : N-1
        for j = 2:M-1 
           Derivatives_u(i,j) = sum(a2(i,:) .* u(:,j)') + sum(b2(j,:) .* u(i,:));
            u(i,j) = u(i,j)+ (B + u(i,j)^2*v(i,j) - (A+1)*u(i,j) + alpha*Derivatives_u(i,j))*dt;
           
            Derivatives_v(i,j) = sum(a2(i,:) .* v(:,j)') + sum(b2(j,:) .* v(i,:));
            v(i,j) = v(i,j) + (A*u(i,j) - u(i,j)^2*v(i,j) + alpha*Derivatives_v(i,j))*dt;
        end
    end
        
         %Update BCs
         [u, v] = updateMatrices(u, v, a1, b1);
     
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
                
    %Storing every value of u
    u_storage(:, :, time) = u;
    v_storage(:, :, time) = v;
    %Time counter
    t = t + dt; 
  end

%Plotting single grid point in time
tt = linspace(0, tau, nsteps);
gridpoint_x = N;
gridpoint_y = M;
values_to_plotu = squeeze(u_storage(gridpoint_x, gridpoint_y, :));
values_to_plotv = squeeze(v_storage(gridpoint_x, gridpoint_y, :));
figure(5)
plot(tt, values_to_plotu);
title(sprintf('U and V in grid-point %d, %d', gridpoint_x, gridpoint_y));
hold on 
plot(tt, values_to_plotv);
%% Save Matrixes
save('u_storage.mat', 'u_storage');
save('v_storage.mat', 'v_storage');

l = matfile('u_storage.mat');
h = matfile('v_storage.mat');
%% Display interesting values
rows = 3;
data = {
    't = 0.5s', u_storage(ceil(N/10), ceil(N*9/10), 0.5/dt), v_storage(ceil(N/10), ceil(N*9/10), 0.5/dt);
    't = 1.0s', u_storage(ceil(N/10), ceil(N*9/10), 1/dt),   v_storage(ceil(N/10), ceil(N*9/10), 1/dt);
    't = 5s', u_storage(ceil(N/10), ceil(N*9/10), 5/dt), v_storage(ceil(N/10), ceil(N*9/10), 5/dt)
       }; 

% Print column headers
fprintf('%10s %10s %10s\n', 'Time', 'u', 'v');

% Print rows of data
for i = 1:rows
   fprintf('%10s %10.2f %10.2f\n', data{i, 1}, data{i, 2}, data{i, 3});
end

%% Create Video
% Initialize video writer
v = VideoWriter('myvideo.avi');
v.FrameRate = 500;  % Set FPS to 30
open(v);

% Get the dimensions of u_storage
dt = 0.002;
v_film = v_storage(:,:,1:2:end);
u_film = u_storage(:,:,1:2:end);
[x, y, time] = size(v_film);

% Create video using surf
figure;
for t = 1:time
    subplot(2,1,1)
    surf(X,Y,u_film(:,:,t));
    titleStr = sprintf('Time = %.2f s', t*dt);  % dt is your time step size
    title(titleStr);
%     axis([0 1 0 1 0 1]);
    subplot(2,1,2)
    surf(X,Y,v_film(:,:,t));
%     axis([0 1 0 1 0 1]);
    titleStr = sprintf('Time = %.2f s', t*dt);  % dt is your time step size
    title(titleStr);
    drawnow;
    
    % Capture the frame
    frame = getframe(gcf);
    
    % Write the frame to the video
    writeVideo(v, frame);
end

% Close the video file
close(v);

