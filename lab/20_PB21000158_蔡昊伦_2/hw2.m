for n = [5,10,20,40]
    n
    i = 0:n;
    x1 = -5 + i/n*10;%第一组数
    x2 = -5*cos((2*i+1)*pi/(2*n+2));%第二组数
    xx = [x1;x2];
    for m = 1:2 
        fprintf("第%d组节点",m);
        x = xx(m,:);%interpolation abscissas
        y = 1./(1+x.*x);%interpolation ordinates
        j = 0:500;
        yi = -5 + j/500*10;%interpolation values 
        s = 0;
        for i = 1:n+1
            z = ones(1, length(yi));
            for k = 1:n+1
                if k ~= i
                    z = z .* (yi - x(k)) / (x(i) - x(k));
                end
            end
            s = s + z * y(i);
        end
        l = s;%lagrange function
        out = max(abs(1./(1+yi.*yi)-l))
    end
end