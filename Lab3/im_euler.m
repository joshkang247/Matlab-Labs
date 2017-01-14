function [y, t] = im_euler(t0,tN,y0,h, f)
    t = t0:h:tN;
    y = zeros(1, length(t));
    y(1) = y0;
    for N = 1:length(t) - 1;
        y(N + 1) = y(N) + 0.5*h*(f(t(N), y(N)) + f(t(N+1), y(N) + f(t(N), y(N))*h));
    end;


    
    