function [u, v] = updateMatrices(u, v, a1, b1)
   [N,M] = size(u);
    
    AXU = a1(1,N)*a1(N,1) - a1(1,1)*a1(N,N);
    AYU = b1(1,M)*b1(M,1) - b1(1,1)*b1(M,M);
    AXV = AXU;
    AYV = AYU;
    
    % For u_{1,j}   SOUTH
    for j = 1:M
        summation = 0;
        for k = 2:N-1
            summation = summation + (a1(N,N) * a1(1,k) - a1(1,N) * a1(N,k)) * u(k,j);
        end
        u(1,j) = summation / AXU;
    end
    
    % For u_{N,j}  NORTH
    for j = 1:M
        summation = 0;
        for k = 2:N-1
            summation = summation + (a1(1,1) * a1(N,k) - a1(N,1) * a1(1,k)) * u(k,j);
        end
        u(N,j) = summation / AXU;
    end
    
        % For u_{i,1}  WEST
    for i = 2:N-1
        summation = 0;
        for k = 2:M-1
            summation = summation + (b1(M,M) * b1(1,k) - b1(1,M) * b1(M,k)) * u(i,k);
        end
        u(i,1) = summation / AYU;
    end

    % For u_{i,M}    EAST
    for i = 2:N-1
        summation = 0;
        for k = 2:M-1
            summation = summation + (b1(1,1) * b1(M,k) - b1(M,1) * b1(1,k)) * u(i,k);
        end
        u(i,M) = summation / AYU;
    end

    % For v_{1,j} 
    for j = 1:M
        summation = 0;
        for k = 2:N-1
            summation = summation + (a1(N,N) * a1(1,k) - a1(1,N) * a1(N,k)) * v(k,j);
        end
        v(1,j) = summation / AXV;
    end

    % For v_{N,j}
    for j = 1:M
        summation = 0;
        for k = 2:N-1
            summation = summation + (a1(1,1) * a1(N,k) - a1(N,1) * a1(1,k)) * v(k,j);
        end
        v(N,j) = summation / AXV;
    end



    % For v_{i,1}
    for i = 2:N-1
        summation = 0;
        for k = 2:M-1
            summation = summation + (b1(M,M) * b1(1,k) - b1(1,M) * b1(M,k)) * v(i,k);
        end
        v(i,1) = summation / AYV;
    end

    % For v_{i,M}
    for i = 2:N-1
        summation = 0;
        for k = 2:M-1
            summation = summation + (b1(1,1) * b1(M,k) - b1(M,1) * b1(1,k)) * v(i,k);
        end
        v(i,M) = summation / AYV;
    end
end