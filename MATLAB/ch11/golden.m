function y = golden(f, rang, delta, n)
a = rang(1);
b = rang(2);
alpha = a + 0.382*(b - a);
beta = a + 0.618*(b - a);

for k = 1:n
    if f(alpha) <= f(beta)
        if beta - a <= delta
            y = alpha;
            return;
        end
        b = beta;
        beta = alpha;
        alpha = a + 0.382*(b - a);
    else
        if b - alpha <= delta
            y = beta;
            return;
        end
        a = alpha;
        alpha = beta;
        beta = a + 0.618*(b - a);
    end
end
end