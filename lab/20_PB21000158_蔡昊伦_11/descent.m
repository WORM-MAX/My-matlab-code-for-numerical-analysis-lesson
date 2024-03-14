function min = descent(f, df, x0, epsilon, n, fileID)
x = x0;
p = df(x);
for k = 1:n
    g = @(y) f(x-y*p);
    [a0, b0] = for_back(g,0,0.01);
    lamda = golden(g,[a0,b0],1e-5,10000);
    x = x - lamda.*p;
    fprintf(fileID,'第%d次迭代：f(x_i)=%.15e, x_1=%.15e,x_2=%.15e\n',k,f(x),x(1),x(2));
    p = df(x);
    if norm(p,2)<epsilon
        min = f(x);
        return;
    end 
end
min = f(x);
end