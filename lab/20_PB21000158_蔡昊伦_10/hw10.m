f = @(x) 0.7*sin(4*pi*x)+sin(10*pi*x);

for n = [128,256]
    x = [0:1/n:1-1/n];
    f_list = f(x);
    g = FFT(f_list);
    fprintf("n = %d\n", n);
    for i = 0:n-1
        fprintf("x_%d=%f, y_%d=%f\n",i, real(g(i+1)), i, imag(g(i+1)));
    end

end

