function [y, x] = rk4(f, a, b, m, y0)
% Euler's method for solving initial value problems
% Input: right-hand side function f(x, y) with two variables, interval [a,
% b], number of grid points m, and initial value y0
% Output: solution y
h = (b - a)/m;
x(1) = a; y(1) = y0;
for i = 1:m
    x(i+1) = x(i) + h;
    k1 = h*f(x(i), y(i));
    k2 = h*f(x(i) + h/2, y(i) + k1/2);
    k3 = h*f(x(i) + h/2, y(i) + k2/2);
    k4 = h*f(x(i+1), y(i) + k3);
    y(i+1) = y(i) +  (k1 + 2*k2 + 2*k3 + k4)/6;
end

% plot
%plot(x, y);

end
