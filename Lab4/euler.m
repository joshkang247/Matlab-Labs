function [x, t] = euler(t0,tN,x0,h, f1, f2)
    t = t0:h:tN;
    x = zeros(2, length(t));
    x(1,1) = x0(1);
    x(2,1) = x0(2);
    
    
    for N = 1:length(t) - 1;
        x(1, N + 1) = x(1, N) + h*(f1(t(N), x(1,N), x(2,N)));
        x(2, N + 1) = x(2, N) + h*(f2(t(N), x(1,N), x(2,N))); 
    end;
    


    
    