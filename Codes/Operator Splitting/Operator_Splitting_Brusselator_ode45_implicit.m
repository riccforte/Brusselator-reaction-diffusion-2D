%% Project CCRE-ATP-APC
clear
clc
close all

%% Grid Generation
%2D Domain
N = 20;    %Points along x direction
M = 20;    %Points along y direction

%Grid Generation
x = linspace(0,1,N);
y = linspace(0,1,M);
dx = x(2)- x(1);
[X, Y] = meshgrid(x, y);

%% Problem Data
%Input Parameters
A = 3.4;
B = 1;              %[-]
alpha = 0.002;      %[-]    

dt = 0.01;         %[s]     %time step
tau = 10;          %[s]     %total time of simulation
nsteps = tau/dt;   %[-]     %total time steps

%% Initial conditions
u = 0.5 + Y;
v = 1 + 5 .* X;

%% Time loop
for t = 1:nsteps
    
    [u, v] = reaction_step(u, v, A, B, dt, N, M);
    [u_new,v_new] = Diffusion(dt,dx,alpha,u,v,N,M);
    u = u_new;
    v = v_new;
    
    % Live plot
    surf(X,Y,u)
    title(['U surface at ' num2str(t*dt)])
    zlabel('V')
    xlabel('X')
    ylabel('Y')
    l = (max(max(u))+0.5);
    axis([0 1 0 1 0 l])
    drawnow    
end
function [u_new,v_new] = Diffusion(dt,dx,alpha,u,v,N,M)

[A] = ConstCoeffImplicitCNMatrix(alpha,dx,N,M,dt);

b = reshape(u(:,:)',[],1);
z = reshape(v(:,:)',[],1);
 
u_new = A\b; % Direct method using LU factorization
v_new = A\z; % Direct method using LU factorization

u_new = reshape(u_new,[N,M])';
v_new = reshape(v_new,[N,M])';

%BCs
u_new(1,:) = u_new(2,:);
u_new(end,:) = u_new(end-1,:);
u_new(:,1) = u_new(:,2);
u_new(:,end) = u_new(:,end-1);
    
v_new(1,:) = v_new(2,:);
v_new(end,:) = v_new(end-1,:);
v_new(:,1) = v_new(:,2);
v_new(:,end) = v_new(:,end-1); 
end
function [u_new, v_new] = reaction_step(u, v, A, B, dt, N, M)
    % Initialize new matrices for u and v
    u_new = zeros(N, M);
    v_new = zeros(N, M);
    
    % Options for ODE solver
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    
    % Reaction step for each grid point
    for i = 1:N
        for j = 1:M
            y0 = [u(i,j); v(i,j)];
            [~, Y] = ode23s(@(t,y) odefun(t, y, A, B), [0, dt], y0, options);
            u_new(i,j) = Y(end, 1);
            v_new(i,j) = Y(end, 2);
        end
    end
end

function dydt = odefun(~, y, A, B)
    u = y(1);
    v = y(2);
    du = B + u^2 * v - (A + 1) * u;
    dv = A * u - u^2 * v;
    dydt = [du; dv];
end

function[A] = ConstCoeffImplicitCNMatrix(k,dxz,i,j,dt)
% This function generates the matrix for solution of 2D diffusion equation
% by Direct implicit method using Crank-Nicholson scheme.

a = 1+(4*k*dt/(2*dxz^2));
c = -k*dt/(2*dxz^2);
ij = i*j;

main_diag = a.*ones(ij,1);
plusi_diag = c.*ones(ij,1);
mini_diag = plusi_diag;

c_unit = c.*ones(i,1);
min1c_unit = c_unit;
min1c_unit(end) = 0;
min1_diag = repmat(min1c_unit,j,1);

plus1c_unit = c_unit;
plus1c_unit(1) = 0;
plus1_diag = repmat(plus1c_unit,j,1);

A = spdiags([mini_diag,min1_diag,main_diag,plus1_diag,plusi_diag],[-i,-1,0,1,i],ij,ij);
spy(A)
end
