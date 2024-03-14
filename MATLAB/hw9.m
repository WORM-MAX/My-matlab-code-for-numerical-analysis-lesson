f = @(x,y) -(x.^2).*y.^2;
y_theory = @(x) 3/(1+x.^3);
y_th = y_theory(1.5);
error = zeros(1,5);
y0 = 3;

fprintf("四阶Runge-Kutta公式的误差和误差阶:\n");
for l=0:4
    h = 0.1/(2^l);
    tspan = 0: h :1.5;
    [~,y] = RungeKutta4(f, tspan, y0);
    error(l+1) = abs(y_th-y(end));
end
for l=1:4
    h = 0.2/(2^l);
    ok = log2(error(l)/error(l+1))/log2(2);
    fprintf("h = %.5f, err = %.15f, ok = %.15f\n", h, error(l), ok);
end

fprintf("\n四阶隐式Adams公式的误差和误差阶:\n");
for l=0:4
    h = 0.1/(2^l);
    tspan = 0: h :1.5;
    [~,y] = Adams4(f, tspan, y0);
    error(l+1) = abs(y_th-y(end));
end
for l=1:4
    h = 0.2/(2^l);
    ok = log2(error(l)/error(l+1))/log2(2);
    fprintf("h = %.5f, err = %.15f, ok = %.15f\n", h, error(l), ok);
end

