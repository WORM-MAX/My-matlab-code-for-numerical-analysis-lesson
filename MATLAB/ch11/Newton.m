function min = Newton(f, df, G, x0, epsilon, n, fileID)
x = x0;
p = -G(x)/df(x);
p = p';
for k = 1:n
    g = @(y) f(x+y*p);
    [a0, b0] = for_back(g,0,1);
    lamda = golden(g,[a0,b0],1e-6,100000);
    x = x + lamda*p;
    fprintf(fileID, '第%d次迭代：f(x_i)=%.15e, x_1=%.15e,x_2=%.15e\n',k,f(x),x(1),x(2));
    p = -df(x)/G(x);
    if norm(p,2)<epsilon
        min = f(x);
        return;
    end 
end
min = f(x);
end