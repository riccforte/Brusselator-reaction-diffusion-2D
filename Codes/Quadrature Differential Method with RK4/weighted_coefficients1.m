function [a] = weighted_coefficients1(N,M,x)
a = zeros(N,M);
PROD=1;
SUM=0;
for i=1:N
    for j=1:N
        if i~=j
        a1=(1/(x(j)-x(i)));
        for k=1:N
            if k~=i && k~=j
              PROD=PROD*(x(i)-x(k))/(x(j)-x(k));    
           end
        end         
        a(i,j)=a1*PROD;
        PROD=1;
        end
        
        if i==j
            for k=1:N
                if k~=i
                   SUM=SUM+(1/(x(i)-x(k))); 
               end
            end
            
            a(i,j)=SUM;
            SUM=0;
        end
    end
end

end