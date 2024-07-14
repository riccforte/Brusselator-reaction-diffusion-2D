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

x = linspace(0,Lx,N);
y = linspace(0,Ly,M);

dx = x(2) - x(1);
dy = y(2) - y(1);

%Creating rectangular Grid
[X ,Y] =meshgrid(x,y);

%% Problem Data
%Input Parameters
A = 3.4;
B = 1;               %[-]
alpha = 0.002;       %[-]     

dt = 0.01;           %[s]     %time step
tau = 10;            %[s]     %total time of simulation
nsteps = tau/dt;     %[-]     %total time steps

%% Initial conditions

u_old = 0.5 + Y;
v_old = 1 + 5 .* X;

%% Implicit transient solver

max_newton_iterations = 1500;
tolerance = 1e-3;

% Start measuring CPU time
tic
for t = dt:dt:tau
    % At each time step, use Newton's method to solve the nonlinear system
    % Assume u_old and v_old are the solutions from the previous time step    
    % Initial guess for the solutions at the new time step
    u = u_old;
    v = v_old;
    
    for newton_iter = 1:max_newton_iterations
        
        % Evaluate the residual of the discretized PDEs
        [R_u, R_v] = evaluate_residuals(u, v, u_old, v_old, A, B, alpha, dx, dy, dt, N, M);
        
        % Inside your time-stepping loop
        J = construct_jacobian(u, v, A, B, alpha, dx, dy, dt, N, M);
    
        % Solve the linear system J*delta = -R to update the solution
        delta = J \ [-R_u(:); -R_v(:)];
        
        % Update the solution vectors
        u = u + reshape(delta(1:numel(u)), size(u));
        v = v + reshape(delta(numel(u)+1:end), size(v));
        
        % Enforce Neumann boundary conditions by setting the flux at the boundaries to zero
        u(1,:) = u(2,:);            % Left boundary
        u(end,:) = u(end-1,:);      % Right boundary
        u(:,1) = u(:,2);            % Bottom boundary
        u(:,end) = u(:,end-1);      % Top boundary

        v(1,:) = v(2,:);            % Left boundary
        v(end,:) = v(end-1,:);      % Right boundary
        v(:,1) = v(:,2);            % Bottom boundary
        v(:,end) = v(:,end-1);      % Top boundary
       
        
        % Check for convergence
        if norm(delta, Inf) < tolerance
              disp(['Converged in ', num2str(newton_iter), ' iterations']);
              
            break;
        end
    end
    
    
    % Check if Newton's method failed to converge
    if newton_iter == max_newton_iterations
        error('Newton''s method did not converge');
    end
    % Update the old solution
    u_old = u;
    v_old = v;
    
    % Plotting
    figure(2)
    mesh(X, Y, v)
    shading interp
    title(['U at time t = ', num2str(t), ' s'])  % Add time indicator to the title
    xlabel('X')
    ylabel('Y')
    zlabel('U')
    drawnow
   

end
% Stop measuring CPU time and display the elapsed time
elapsed_time = toc;
disp(['Elapsed CPU time: ', num2str(elapsed_time), ' seconds']);
function [R_u, R_v] = evaluate_residuals(u, v, u_old, v_old, A, B, alpha, dx, dy, dt, N, M)
    % Initialize the residual arrays
    R_u = zeros(N, M);
    R_v = zeros(N, M);

    % Interior points
    for i = 2:N-1
        for j = 2:M-1
            % Discretize the Laplacian for u and v using central differences
            lap_u = (u(i+1,j) - 2*u(i,j) + u(i-1,j)) / dx^2 + ...
                    (u(i,j+1) - 2*u(i,j) + u(i,j-1)) / dy^2;
            lap_v = (v(i+1,j) - 2*v(i,j) + v(i-1,j)) / dx^2 + ...
                    (v(i,j+1) - 2*v(i,j) + v(i,j-1)) / dy^2;
            
            % Reaction terms
            react_u = B + u(i,j)^2 * v(i,j) - (A + 1) * u(i,j);
            react_v = A * u(i,j) - u(i,j)^2 * v(i,j);
            
            % Combine reaction and diffusion for the time derivative
            R_u(i,j) = (u(i,j) - u_old(i,j)) / dt - react_u - alpha * lap_u;
            R_v(i,j) = (v(i,j) - v_old(i,j)) / dt - react_v - alpha * lap_v;
        end
    end
    
    % Neumann boundary conditions (zero flux)
    R_u(1,:) = u(2,:) - u(1,:);       % Left boundary
    R_u(N,:) = u(N,:) - u(N-1,:);     % Right boundary
    R_u(:,1) = u(:,2) - u(:,1);       % Bottom boundary
    R_u(:,M) = u(:,M) - u(:,M-1);     % Top boundary

    R_v(1,:) = v(2,:) - v(1,:);       % Left boundary
    R_v(N,:) = v(N,:) - v(N-1,:);     % Right boundary
    R_v(:,1) = v(:,2) - v(:,1);       % Bottom boundary
    R_v(:,M) = v(:,M) - v(:,M-1);     % Top boundary

    % Flatten the residuals to vectors
    R_u = R_u(:);
    R_v = R_v(:);
end

function J = construct_jacobian(u, v, A, B, alpha, dx, dy, dt, N, M)
    % Initialize the Jacobian matrix as a sparse matrix
    % There are N*M equations for u and N*M equations for v
    J = spalloc(2*N*M, 2*N*M, 5*2*N*M);

    % Helper function to convert (i, j) index to linear index for u and v
    idx_u = @(i, j) (j-1)*N + i;
    idx_v = @(i, j) N*M + (j-1)*N + i;

    % Compute the Jacobian matrix
    for i = 1:N
        for j = 1:M
            % Linear indices for (i,j)
            ind_u = idx_u(i, j);
            ind_v = idx_v(i, j);

            % Diagonal entries for u and v
            J(ind_u, ind_u) = 1/dt - (A + 1) - 2*u(i,j)*v(i,j) + 2*alpha*(1/dx^2 + 1/dy^2);
            J(ind_v, ind_v) = 1/dt + 2*u(i,j)*v(i,j) + 2*alpha*(1/dx^2 + 1/dy^2);

            % Off-diagonal entries for u's equation
            J(ind_u, ind_v) = u(i,j)^2;

            % Off-diagonal entries for v's equation
            J(ind_v, ind_u) = -A + 2*u(i,j)*v(i,j);

            % Now, handle the Laplacian operator terms for u
            if i > 1
                J(ind_u, idx_u(i-1, j)) = -alpha/dx^2;
            end
            if i < N
                J(ind_u, idx_u(i+1, j)) = -alpha/dx^2;
            end
            if j > 1
                J(ind_u, idx_u(i, j-1)) = -alpha/dy^2;
            end
            if j < M
                J(ind_u, idx_u(i, j+1)) = -alpha/dy^2;
            end

            % Handle the Laplacian operator terms for v
            if i > 1
                J(ind_v, idx_v(i-1, j)) = -alpha/dx^2;
            end
            if i < N
                J(ind_v, idx_v(i+1, j)) = -alpha/dx^2;
            end
            if j > 1
                J(ind_v, idx_v(i, j-1)) = -alpha/dy^2;
            end
            if j < M
                J(ind_v, idx_v(i, j+1)) = -alpha/dy^2;
            end
        end
    end

    % Correct the Jacobian entries at the boundaries if Neumann conditions apply
    for i = 1:N
        for j = 1:M
            if i == 1 || i == N || j == 1 || j == M
                % Modify the Jacobian for Neumann boundary conditions
                % For example, if the derivative at the boundary should be zero
                % then the difference between the boundary and the adjacent interior
                % point should also be zero. This means we add an additional term
                % to the diagonal of the boundary points equal to what would have
                % been subtracted for the missing neighbor.
                % If i == 1, there is no left neighbor
                J(idx_u(i,j), idx_u(i,j)) = J(idx_u(i,j), idx_u(i,j)) + (i == 1)*alpha/dx^2;
                % If i == N, there is no right neighbor
                J(idx_u(i,j), idx_u(i,j)) = J(idx_u(i,j), idx_u(i,j)) + (i == N)*alpha/dx^2;
                % If j == 1, there is no neighbor below
                J(idx_u(i,j), idx_u(i,j)) = J(idx_u(i,j), idx_u(i,j)) + (j == 1)*alpha/dy^2;
                % If j == M, there is no neighbor above
                J(idx_u(i,j), idx_u(i,j)) = J(idx_u(i,j), idx_u(i,j)) + (j == M)*alpha/dy^2;
                
                % Repeat for v
                J(idx_v(i,j), idx_v(i,j)) = J(idx_v(i,j), idx_v(i,j)) + (i == 1)*alpha/dx^2;
                J(idx_v(i,j), idx_v(i,j)) = J(idx_v(i,j), idx_v(i,j)) + (i == N)*alpha/dx^2;
                J(idx_v(i,j), idx_v(i,j)) = J(idx_v(i,j), idx_v(i,j)) + (j == 1)*alpha/dy^2;
                J(idx_v(i,j), idx_v(i,j)) = J(idx_v(i,j), idx_v(i,j)) + (j == M)*alpha/dy^2;
            end
        end
    end
end
