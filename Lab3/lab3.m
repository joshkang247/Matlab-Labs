%% ODE Lab: Creating your own ODE solver in MATLAB
%
% In this lab, you will write your own ODE solver for the Heun method (or 
% Improved Euler method), and compare its results to those of |ode45|.
%
% You will also learn how to write a function in a separate m-file and 
% execute it.
% 
% Opening the m-file lab3.m in the MATLAB editor, step through each
% part using cell mode to see the results.  Compare the output with the
% PDF, which was generated from this m-file.
%
% There are six (6) exercises in this lab that are to be handed in at the
% end of the lab.  Write your solutions in the template, including
% appropriate descriptions in each step. Save the .m files and submit them 
% online using Blackboard.
%
% MAT292, Fall 2016, Sinnamon & Sousa, modified from
% MAT292, Fall 2011, Hart & Pym

%% Student Information
%
% Student Name: Jiaxi Kang
%
% Student Number: 1002413328
%

%% Creating new functions using m-files.
%  
% Create a new function in a separate m-file:
%
% Specifics:  Create a text file with the file name f.m
% with the following lines of code (text):
%
%  function y = f(a,b,c) 
%  y = a+b+c;
%
% Now matlab can call the new function f (which simply accepts 3 numbers
% and adds them together).  
% To see how this works, type the following in the matlab command window:
% sum = f(1,2,3)

%% Exercise 1
%
% Objective: Write your own ODE solver (using the Heun/Improved Euler
% Method).
%
% Details: This m-file should be a function which accepts as variables 
% (t0,tN,y0,h), where t0 and tN are the start and end points of the 
% interval on which to solve the ODE, y0 is the initial condition of the
% ODE, and h is the stepsize.  You may also want to pass the function into
% the ODE the way |ode45| does (check MATLAB lab 2).
%
% Note: you will need to use a loop to do this exercise.  
% You will also need to recall the Heun/Improved Euler algorithm learned in lectures.  


%

%% Exercise 2
%
% Objective: Compare Heun with |ode45|.
%
% Specifics:  For the initial-value problems of exercises 1, 3, 4 and 5 
% from MATLAB lab 2, approximate the solutions with your function from
% exercise 1 (Improved Euler Method).
% Plot the grahs of your Improved Euler Approximation with the |ODE45| 
% approximation.
%
% Comment on any major differences, or the lack thereof. You do not need
% to reproduce all the code here. Simply make note of any differences of
% note exercise by exercise.


%exercise 1
% f = @(t, y) y.*tan(t) + sin(t);
% t0 = 0;
% y0 = -0.5;
% tN = pi;
% 
% [y, t] = im_euler(t0, tN, y0, 0.1, f);
% soln = ode45(f, [t0, tN], y0);
% plot(t, y, soln.x, soln.y);
% xlabel('t');
% legend('Improved euler', 'ode45');

%as h is very small, improved euler is almost identical to ode45. However,
%as h is slightly large it appears there is a kink pi/2. This is caused by
%tan(pi/2) = infinity. 

%exercise 3
% yexact1 = -cos(soln.x) ./ 2;
% yexact2 = -cos(t) ./ 2;
% err1 = abs(yexact1 - soln.y);
% err2 = abs(yexact2 - y);
% fprintf('maximum error of ode: %g\n', max(err1));
% fprintf('maximum error of improved euler: %g\n', max(err2));
% 
% semilogy(t, err2);
% hold on;
% semilogy(soln.x, err1);
% xlabel('t');
% ylabel('error');

%shape of error is similar, both peaking at pi/2. However, as h is very
%small, the values of the error are similar for both ode45 and improved
%euler. However, as h becomes larger, errors for improved euler are much
%higher.


%exercise 4
% f = @(t, y) 1 ./  (y.^2);
% t0 = 1;
% y0 = 1;
% tN = 10;
% 
% [y, t] = im_euler(t0, tN, y0, 0.01, f);
% soln = ode45(f, [t0, tN], y0);
% plot(t, y, soln.x, soln.y);
% xlabel('t');
% legend('Improved euler', 'ode45');
% 
% yexact1 = (-2 + 3.*(soln.x)).^(1/3);
% yexact2 = (-2 + 3.*(t)).^(1/3);
% err1 = abs(yexact1 - soln.y);
% err2 = abs(yexact2 - y);
% fprintf('maximum error of ode: %g\n', max(err1));
% fprintf('maximum error of improved euler: %g\n', max(err2));

%in this exercise a high h such as 1 will produce greater errors than the
%ode. However, an h such as 0.01 will produce a much smaller error than the
%ode as seen from err1 and err2.


%exercise 5
% f = @(t, y) 1 - t.*y ./ 2;
% t0 = 0;
% y0 = -1;
% tN = 10;
% 
% [y, t] = im_euler(t0, tN, y0, 1, f);
% soln = ode45(f, [t0, tN], y0);
% plot(t, y, soln.x, soln.y);
% xlabel('t');
% legend('Improved euler', 'ode45');

%ode45 shows a curve with many kinks, whereas improved euler with a low
%h=0.01, shows a smooth curve. However, if h = 1, we get a very different
%curve than the solution from ode45. In this case ode45 solution at t=10 is
%approximately 0.2 whereas improved euler at h=1 at t=10 is almost -230.

%% Exercise 3
%
% Objective: Use Euler method and verify that the estimate for the global
% truncation error studied in class is valid.
%
% Details: 
%
% (a) Use Euler method (you can use
% euler.m from iode) to solve the IVP
%
%  |y' = 2 t sqrt( 1 - y^2 )  ,  y(0) = 1/sqrt(3)|
%
% from |t=0| to |t=0.5|.

f = @(t, y) 2.*t.*sqrt(1-y.^2);
t0 = 0;
tN = 0.5;
h = 0.1;
y0 = 1/sqrt(3);
t = t0:h:tN;
y = euler(f, y0, t);
y(length(t)) %0.7587

% (b) Calculate the solution of the IVP and evaluate it at |t=0.5|.

exact = @(t) sin(t.^2 + asin(1/sqrt(3)));
exact(0.5) %0.7614

% (c) In lecture, we learned how to estimate the global truncation 
% error for the Euler method. Write the estimate (as a comment), specifying
% the meaning of all the constants used.
% En = (kh/2M)*(e^(Mnh) -1)
% k = A + MB
% A = max(df/dt)
% B = max(f)
% C = max(df/dy)
%h = step-size


% (d) Compute the error estimate for |t=0.5| and compare with the actual
% error.

% f(0.5, y) = sqrt(1-y^2), max = 1
% df/dt(0.5, y) = 2*sqrt(1-y^2) max = 2
% df/dy(0.5, y) = (-y(1-y^2)^-0.5) max = 1.1745
M = 1.1745;
A = 2;
B = 1;
K = A + M*B;
n = 0.5/h;
err = ((K*h)/(2*M))*(exp(M*n*h) - 1)

err_exact = exact(0.5) - y(length(t))

%actual error is much smaller than estimated

% (e) Change the time step and compare the new error estimate with the
% actual error. Comment on how it confirms the global truncation error
% studied in lectures.

% it is known that Eulers method error is proportional to h to the power of
% 1. However, the actual error appears also to be proportional to h. As h
% changes from 0.1 to 0.01 actual error changes from 0.0295 to 0.0027,
% approximately a linear relationship. Thus we can confirm that Euler's
% method is actually of order one. 


%% Adaptive Step Size
%
% As mentioned in MATLAB assignment 2, the step size in |ode45| is adapted to a
% specific error tolerance.
%
% The idea of adaptive step size is to change the step size |h| to a
% smaller number whenever the derivative of the solution changes quickly.
% This is done by evaluating f(t,y) and checking how it changes from one
% iteration to the next.

%% Exercise 4
%
% Objective: Create an Adaptive Euler method, with an adaptive step size |h|.
%
% Details: Create an m-file which accepts the variables |(t0,tN,y0,h)|, as 
% in exercise 1, where h is the maximum step size. You may also want to 
% pass the function into the ODE the way |ode45| does.
%
% Use Euler's method from |iode| and change it to include the following:
%
% (a) Approximation of the derivative of |f(t,y(t))|:
%
%  |d/dt ( f(t,y(t) ) = [ df/dt + df/dy . f ] (t,y(t))|
%
% by 
%
%  |{ [ f(t[n+1],y[n]) - f(t[n],y[n]) ] + [ f(t[n],y[n+1]) - f(t[n],y[n]) ] } / h|
%  Explain why we use this formula.

% d/dt(f(t, y(t)) can be approximated by a change in f in a small amount of
% time, h. The first part of the equation (f(t[n+1],y[n]) - f(t[n],y[n]))/h
% determines the change in f due to a change in time, which approximates
% df/dt. The second part of the equation (f(t[n],y[n+1]) - f(t[n],y[n]))/h
% represents df/dy. This equation can then be used numerically in matlab.
% By knowing the second derivative of y, it will determine how fast the
% slow is changing and how much we need to adjust h to reduce the error.
% With a high second derivative, a small h is needed to gain an accurate
% approximation. On the other hand, with a low second derivative, a larger
% size h is sufficient for an accurate approximation. 
%
%
% (b) At every iteration compute this approximation and check whether or
% not it is greater than 1.1 times that of the previous iteration.

% (c) If it is, then re-compute that iteration with |h/2| instead of |h|
% (decrease the time step by half).
%
% When you program this, make sure to include a way to avoid infinite
% loops.
%
% (d) Continue until the approximation reaches the predetermined |tN|.

%% Exercise 5
%
% Objective: Compare Euler to your Adaptive Euler method.
%
% Details: Consider the IVP from exercise 3.
%
% (a) Use Euler method to approximate the solution from |t=0| to |t=0.75|
% with |h=0.025|.

f = @(t, y) 2.*t.*sqrt(1-y.^2);
t0 = 0;
tN = 0.75;
h = 0.025;
y0 = 1/sqrt(3);

t1 = t0:h:tN;
y1 = euler(f, y0, t1);
% (b) Use your Adaptive Euler method to approximate the solution from |t=0| 
% to |t=0.75| with initial |h=0.025|.

[y2, t2] = ad_euler(t0,tN,y0,h,f);


% (c) Plot both approximations together with the exact solution.
exact = @(t) sin(t.^2 + asin(1/sqrt(3)));
plot(t1, y1, t2, y2, t2, exact(t2));
legend('Euler', 'Adaptive Euler', 'Exact');


%% Exercise 6
%
% Objective: Problems with Numerical Methods.
%
% Details: Consider the IVP from Exercise 3 (and 5).
% 
% (a) From the two approximations calculated in exercise 5, which one is
% closer to the actual solution? Explain why.

err_adaptive = abs(exact(t2) - y2);
err_euler = abs(exact(t1) - y1);

err_a = sum(err_adaptive)/length(err_adaptive); %0.0049
err_e = sum(err_euler)/length(err_euler); %0.0047

%To determine which is is more accurate, it is necessary to calculate the
%error at each point for each method, and sum them up. However, the number
%of points for each method are different, so it is necessary to normalize
%them by the number of steps required for each method. Surprisinly, it
%turns out that the adaptive method is less accurate than the euler method
%by a very small amount. However, as we check the error when tN is about 1
%instead of 0.75, we notice that euler is actually less accurate than
%adaptive euler. This can be caused by the choice of the ratio to half the
%time step by 2. In this case, the ratio of 1.1 may not be appropriate for
%t(0, 0.75).



% (b) Plot the Exact Solution (done in 3.b), the Euler's Approximation
% (done in 3.a) and the Adaptive Euler's Approximation (done in 5) from
% |t=0| to |t=1.5|.

f = @(t, y) 2.*t.*sqrt(1-y.^2);
t0 = 0;
tN = 1.5;
h = 0.025;
y0 = 1/sqrt(3);

t1 = t0:h:tN;
y1 = euler(f, y0, t1);

[y2, t2] = ad_euler(t0,tN,y0,h,f);
exact = @(t) sin(t.^2 + asin(1/sqrt(3)));
plot(t1, y1, t2, y2, t1, exact(t1));
legend('Euler', 'Adaptive Euler', 'Exact');

% (c) Notice how the exact solution and the approximations become very
% different. Why is that? Write your answer as a comment.

% The cause of the difference begins at t = 1, y = 1. Although the exact
% solution is continous, the derivative is not. At y = 1, the derivative
% becomes imaginery and the plot only shows the real part of the solution.
% This can be shown by simply evaluating the functions calculated by both
% methods at a t >= 1;
