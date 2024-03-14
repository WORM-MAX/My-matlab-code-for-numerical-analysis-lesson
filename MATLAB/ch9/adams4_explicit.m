function [y, x] = adams4_explicit(f, a, b, m, y0)
% Adams 4-order explicit method for solving initial value problems
% Input: right-hand side  function(x, y) with two variables, interval [a,
% b], number of grid points m, and initial value y0
% Output: solution y
h = (b - a)/m;
x = [a:h:b];
[y_initial,] = rk4(f, a, b, m, y0);
y = y_initial;
for i = 4:m
    y(i+1) = y(i) +  h/24*[55*f(x(i), y(i))-59*f(x(i-1), y(i-1))+37*f(x(i-2), y(i-2))-9*f(x(i-3), y(i-3))];
end

% plot
%plot(x, y);

end
