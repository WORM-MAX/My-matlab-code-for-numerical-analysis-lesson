function [y, x] = adams4_implicit(f, a, b, m, y0)
% Adams 4-order implicit method for solving initial value problems
% Input: right-hand side function f(x, y) with two variables, interval [a,
% b], number of grid points m, and initial value y0
% Output: solution y
h = (b - a)/m;
x = [a:h:b];

%四阶RK起步
[y_initial,] = rk4(f, a, b, m, y0);
y = y_initial;

%起步校正
y(4) = y(3) + h/24*[9*f(x(4), y(4))+19*f(x(3), y(3))-5*f(x(2), y(2))+f(x(1), y(1))];

for i = 4:m
    y(i+1) = y(i) +  h/24*[55*f(x(i), y(i))-59*f(x(i-1), y(i-1))+37*f(x(i-2), y(i-2))-9*f(x(i-3), y(i-3))];%预估
    y(i+1) = y(i) + h/24*[9*f(x(i+1), y(i+1))+19*f(x(i), y(i))-5*f(x(i-1), y(i-1))+f(x(i-2), y(i-2))];%校正
end

% plot
%plot(x, y);

end
