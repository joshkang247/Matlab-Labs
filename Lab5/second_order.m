function [y, t, v] = second_order(t0,tN,y0, y1,h, f)
    t = t0:h:tN;
    
    y = zeros(1,length(t));
    v = zeros(1,length(t));
    v(1) = y1;
    
    y(1) = y0;
    y(2) = y0 + h*y1;
    
    
    
    for n = 2:length(t) - 1;
        v(n) = (y(n) - y(n-1))/h;
        y(n + 1) = (h^2)*f(t(n),y(n),v(n)) + 2*y(n) - y(n-1);
    end;
    
    v(length(t)) = (y(length(t)) - y(length(t)-1))/h;
    