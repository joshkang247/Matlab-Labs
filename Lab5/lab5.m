%% Second-Order Lab: Second-Order Linear DEs in MATLAB
%
%
% In this lab, you will learn how to use |iode| to plot solutions of 
% second-order ODEs. You will also learn to classify the behaviour of 
% different types of solutions.
%
% Moreover, you will write your own Second-Order ODE system solver, and 
% compare its results to those of |iode|.
%
% Opening the m-file lab5.m in the MATLAB editor, step through each
% part using cell mode to see the results.  Compare the output with the
% PDF, which was generated from this m-file.
%
% There are seven (7) exercises in this lab that are to be handed in at the
% end of the lab.  Write your solutions in the template, including
% appropriate descriptions in each step. Save the m-files and submit them 
% online using Blackboard.
%
% MAT292, Fall 2016, Sinnamon & Sousa, based on
% Spring 2010, Jerrard


%% Student Information
%
%  Student Name: Jiaxi Kang
%
%  Student Number: 1002413328
%

%% Iode for Second-Order Linear DEs with constant coefficients
% 
% In the |iode| menu, select the |Second order linear ODEs| module.
% It opens with a default DE and a default forcing function |f(t) =
% cos(2t)|. The forcing function can be plotted along with the solution by
% choosing |Show forcing function| from the |Options| menu.
%
% Use this module to easily plot solutions to this kind of equations.
%
% There are three methods to input the initial conditions:
%
% Method 1. Enter the values for |t0|, |x(t0)|, and |x'(t0)| into the
% |Initial conditions| boxes, and then click |Plot solution|.
%
% Method 2. Enter the desired slope |x'(t0)| into the appropriate into the
% |Initial conditions| box, and then click on the graph at the point
% |(t0,x(t0))| where you want the solution to start.
%
% Method 3. Press down the left mouse button at the desired point
% |(t0,x(t0))| and drag the mouse a short distance at the desired slope
% |x'(t0)|. When you release the mouse button, |iode| will plot the
% solution.

%% Growth and Decay Concepts
%
% We want to classify different kinds of behaviour of the solutions. We say
% that a solution:
%
% |grows| if its magnitude tends to |infinity| for large values of |t|, that
% is, if either the solution tends to $+\infty$ or $-\infty$,
%
% |decays| if its magnitude converges to |0| for large values of |t|,
%
% |decays while oscillating| if it keeps changing sign for large values of
% |t| and the amplitude of the oscillation tends to zero,
%
% |grows while oscillating| if it keeps changing sign for large values of
% |t| and the amplitude of the oscillation tends to infinity.

%% Example
%
t = 0:0.1:10;

% Example 1
figure();
y1 = exp(t);
plot(t,y1)

% Annotate the figure
xlabel('t');
ylabel('f1(t)');
title('The function e^t grows');
legend('f1(t)=e^t');



% Example 2
figure();
y2 = -exp(t);
plot(t,y2)

% Annotate the figure
xlabel('t');
ylabel('f2(t)');
title('The function e^{-t} grows');
legend('f2(t)=-e^{t}');



% Example 3
figure();
y3 = exp(-t);
plot(t,y3)

% Annotate the figure
xlabel('t');
ylabel('f3(t)');
title('The function e^{-t} decays');
legend('f3(t)=e^{-t}');



% Example 4
figure();
y4 = exp(-t).*cos(t);
plot(t,y4)

% Annotate the figure
xlabel('t');
ylabel('f4(t)');
title('The function e^{-t}*cos(t) decays while oscillating');
legend('f4(t)=e^{-t}*cos(t)');



% Example 5
figure();
y5 = exp(t).*sin(2*t);
plot(t,y5)

% Annotate the figure
xlabel('t');
ylabel('f5(t)');
title('The function e^{t}sin(2t) grows while oscillating');
legend('f5(t)=e^{t}sin(2t)');



% Example 6
figure();
y6 = sin(3*t);
plot(t,y6)

% Annotate the figure
xlabel('t');
ylabel('f6(t)');
title('The function sin(3t) neither decays nor grows, it just oscillates');
legend('f6(t)=sin(3t)');



% |Remark.| A function which |grows while oscillating| doesn't |grow|,
% because it keeps changing sign, so it neither tends to $+\infty$ nor to
% $-\infty$.

%% Exercise 1
%
% Objective: Use |iode| to solve second-order linear DEs. And classify
% them.
%
% Details: Consider the ODE:
%
%  4y'' + 4 y' + 17 y = 0
%
% (a) Use |iode| to plot six (6) numerical solutions of this equation with
% "random" initial data (use Method 3 above) and press-and-drag at various
% initial points, with some of the slopes being positive and some negative)
%
% Use only initial points in the part of the window where |0<t<1| and
% |-1<x<1| and take all initial slopes between |-3| and |+3|.
%
% Change the window to [0,10]x[-3,3].
%
% Annotate each solution on the plot title and write whether it decays, grows, or
% decays/grows while oscillating.
%
% (b) Based on the results of (a), state what percentage of solutions 
% decay, grow, grow while oscillating, or decay while oscillating.
%
% (c) Solve the DE and write the exact solution. Explain why this justifies
% your answer in (b).

%% Exercise 2
%
% Consider the ODE:
%
%  y'' + sqrt(3) y' - y/4 = 0
%
% Repeat (a) and (b) from Exercise 1 with this DE.

%% Exercise 3
%
% Consider the ODE:
%
%  y'' + sqrt(3) y' + y/4 = 0
%
% Repeat (a) and (b) from Exercise 1 with this DE.


%% Example
%
% Consider the ODE:
%
%  y'' + 2y' + 10y = 0
% 
% The solution is
%
%  y(t) = e^(-t) ( c1 cos(3t) + c2 sin(3t) )
%
% From this, it is easy to see that all solutions decay while oscillating.
%
% Similarly, for the equation 
%
%  y'' - 2y' + 10y = 0
% 
% The solution is
%
%  y(t) = e^t ( c3 cos(3t) + c4 sin(3t) )
%
% which grows while oscillating.

%% Exercise 4
%
% Consider the fourth-order ODE:
%
%  y'''' + 2 y''' + 6 y'' + 2 y' + 5 y = 0
%
% (a) Find the general solution for this problem (you can use MATLAB to
% find the roots of the characteristic equation).
%
%
% (b) Predict what percentage of solutions with random initial data will 
% grow, decay, grow while oscillating, and decay while oscillating. Explain.

%% Exercise 5
%
% Objective: Classify equations given the roots of the characteristic
% equation.
%
% Details: Your answer can consist of just a short sentence, as |grows| or 
% |decays while oscillating|.
%
% Consider a second-order linear constant coefficient homogeneous DE with
% |r1| and |r2| as roots of the characteristic equation.
%
% Summarize your conclusions about the behaviour of solutions for randomly
% chosen initial data when.
%
% (a) |0 < r1 < r2|
% Unstable node: 100% of solutions will grow.
%
% (b) |r1 < 0 < r2|
% Saddle point: 99% of solutions will grow, 1% will decay
%
% (c) |r1 < r2 < 0|
% Stable node: 100% of solutions will decay
% 
% (d) |r1 = alpha + beta i| and  |r2 = alpha - beta i| and |alpha < 0|
% Spiral sink: 100% of solutions will decay while oscillating 
%
% (e) |r1 = alpha + beta i| and  |r2 = alpha - beta i| and |alpha = 0|
% Center: 100% of solutions will neither grow nor decay, will only
% oscillate
% 
% (f) |r1 = alpha + beta i| and  |r2 = alpha - beta i| and |alpha > 0|
% Spiral source: 100% of solutions will grow while oscillating

%% Numerical Methods for Second-Order ODEs
%
% One way to create a numerical method for second-order ODEs is to
% approximate derivatives with finite differences in the same way of the
% Euler method.
%
% This means that we approximate the first derivative by:
%
%  y'(t[n]) ~ ( y[n] - y[n-1] ) / h
%
% and
%   
%  y''(t[n])  ~ ( y'(t[n+1]) - y'(t[n]) ) / h ~ ( y[n+1] - 2 y[n] + y[n-1]
%  ) / ( h^2 )
%
% By writing these approximations into the ODE, we obtain a method to get
% |y[n+1]| from the previous two steps |y[n]| and |y[n-1]|.
%
% The method for approximating solutions is:
%
% 1. Start with |y[0]=y0|
%
% 2. Then we need to get |y[1]|, but we can't use the method, because we
% don't have two iterations |y[0]| and |y[-1]|(!!). So we use Euler to get
%
%  y[1] = y0 + y1 h
%
% 3. Use the method described above to get |y[n]| for |n=2,3,...|.

%% Exercise 6
%
% Objective: Write your own second-order ODE solver.
%
% Details: Consider the second-order ODE
%
%  y'' + p(t) y' + q(t) y = g(t)
% 
% Write a second-order ODE solver using the method described above.
% 
% This m-file should be a function which accepts as variables 
% (t0,tN,y0, y1,h), where t0 and tN are the start and end points of the 
% interval on which to solve the ODE, y0, y1 are the initial conditions of
% the ODE, and h is the stepsize.  You may also want to pass the functions 
% into the ODE the way |ode45| does (check MATLAB lab 2).
%
% Note: you will need to use a loop to do this exercise.  

%% Exercise 7
%
% Objective: Compare your method with |iode|
%
% Details: Use |iode| to plot the solution of the ODE in Exercise 2 with
% the initial conditions 
%
%  y(0) = 1
%  y'(0) = (sqrt(3)-2)/2
%
% Without removing the figure window, plot your solution (in a different colour), 
% which will be plotted in the same graph.
%
% Comment on any major differences, or the lack thereof. You do not need
% to reproduce all the code here. Simply make note of any differences of
% note.
h = 0.01;
y0 = 1;
t0 = 0;
tN = 10;
y1 = ((sqrt(3)-2)/2);
%y'' + sqrt(3) y' - y/4 = 0
f = @(t, y, v) (1/4).*y - sqrt(3).*v;
[y,t,v] = second_order(t0,tN,y0, y1,h, f);

plot(t, y);
xlabel('t');
ylabel('f(t)');

%We notice that iode and the second_order solver overlap almost exactly
%when seen on the plot on iode (refer to pdf). For both solvers, the
%solution grows.




