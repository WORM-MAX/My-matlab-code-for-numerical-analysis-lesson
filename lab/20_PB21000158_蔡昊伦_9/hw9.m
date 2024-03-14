f = @(x,y) -x^2*y^2;
a = 0;
b = 1.5;
y0 = 3;
y_standard = @(x) 3/(1+x^3);

error_1 = zeros(1,4);
error_2 = zeros(1,4);
for l = 0:3
    m = 15*2^l;
    [y,] = rk4(f, a, b, m, y0);
    error_1(l+1) = abs(y(m+1) - y_standard(1.5));
end

for l = 0:3
    m = 15*2^l;
    [y,] = adams4_implicit(f, a, b, m, y0);
    error_2(l+1) = abs(y(m+1) - y_standard(1.5));
end

%误差阶计算
order_1=zeros(1,4);
order_2=zeros(1,4);
for i=2:4
    order_1(i)=log(error_1(i-1)/error_1(i))/log(2);
    order_2(i)=log(error_2(i-1)/error_2(i))/log(2);
end

%打印
fprintf("四阶Runge-Kutta公式的误差和误差阶：\n");
for i=1:4
    fprintf("h = %.5f, err = %.15f, ok = %.15f\n", 0.1/2^(i-1), error_1(i), order_1(i));
end    

fprintf("四阶阶隐式Adams公式的误差和误差阶：\n");
for i=1:4
    fprintf("h = %.5f, err = %.15f, ok = %.15f\n", 0.1/2^(i-1), error_2(i), order_2(i));
end    