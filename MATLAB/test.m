function yy = f(xx)
yy = 1./(1+xx.*xx)
end

for n = [5,10,20,40]
    i = 1:n
    x = -5 + i/n*10
    y = f(x)
end