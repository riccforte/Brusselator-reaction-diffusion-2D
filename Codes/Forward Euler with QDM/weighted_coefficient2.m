function a2 = weighted_coefficient2(N, M, a1, x)

a2o = zeros(N, M);
difference = 1;
threshold = 1e-8;
while difference > threshold
    a2 = a2o;
    for i = 1 : N
        for j = 1: N
            if i ~= j
                a2(i,j) = 2 *( a1(i,j) * a1(i,i) - a1(i,j) / (x(i) - x(j)));
            else           
               a2(i,i) = -sum(a2o(i, [1:i-1, i+1:end]));  % Summing excluding the diagonal element
            end
        end
    end
    difference = norm(a2 - a2o, 'inf');
    a2o = a2;
  
end

end