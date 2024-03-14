function y = fibonacci(f, a0, b0, delta, n)
F = @(k) (((1+sqrt(5))/2)^(k+1) - ((1-sqrt(5))/2)^(k+1));

a = a0;b = b0;

alpha = a + (1-F(n)/F(n+1))*(b - a);
beta = a + (F(n)/F(n+1))*(b - a);

for k = 1:n
    alpha = a + (1-F(n-k)/F(n-k+1))*(b - a);
    beta = a + (F(n-k)/F(n-k+1))*(b - a);
    
    if f(alpha) <= f(beta)
        b = beta;
    else
        a = alpha;
    end
    a
    b
    if b - a <= delta
        y = a;
        return;
    end
end

end