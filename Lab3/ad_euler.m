function [y, t] = ad_euler(t0,tN,y0,h_init,f)
    h = h_init;
    
    t = [t0];
    
    y = [y0];
    
    y_new = y(length(y)) + h*f(t(length(t)), y(length(y)));
    t_new = t0 + h;
    y = [y y_new];
    t = [t t_new];
    g1 = (f(t(2), y(1)) - f(t(1), y(1)) + f(t(1), y(2)) - f(t(1), y(1)))/h;
    
   
    while t(length(t)) < tN;
        h = h_init;
        y_new = y(length(y)) + h*(f(t(length(t)), y(length(y))));
        t_new = t(length(t)) + h;
        
        
        
        g2 = (f((t_new), y(length(y))) - f(t(length(t)), y(length(y))) + f(t(length(t)), y_new) - f(t(length(t)), y(length(y))))/h;
        while (abs(g2/g1) > 1.1) && (g1 ~= 0);
            if h < h_init/2^4;
                break;
            end
            h = h/2;
            y_new = y(length(y)) + h*(f(t(length(t)), y(length(y))));
            t_new = t(length(t)) + h;
            g2 = (f((t_new), y(length(y))) - f(t(length(t)), y(length(y))) + f(t(length(t)), y_new) - f(t(length(t)), y(length(y))))/h;          
        end
        
        g1 = g2;
        y = [y y_new];
        t = [t t_new];
    end
    

    
    
