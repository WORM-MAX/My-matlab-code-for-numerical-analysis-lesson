f = @(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2;
df = @(x) [400*x(1)*(x(1)^2-x(2))+2*(x(1)-1), 200*(x(2)-x(1)^2)];
G = @(x) [2+100*(8*x(1)^2-4*(-x(1)^2+x(2))),-400*x(1)%hessian matrix of f
    -400*x(1),200];
fileID = fopen('data.txt','w');
fprintf(fileID,'最速下降法:\n');
y = descent(f, df, [0 0], 1e-4, 10000,fileID);
fprintf(fileID, "牛顿法\n");
z = Newton(f, df, G, [0 0], 1e-4, 20, fileID);