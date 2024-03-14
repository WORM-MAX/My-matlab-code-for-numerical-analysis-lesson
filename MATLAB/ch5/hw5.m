%function and derivative
f = @(x) x^3/3-x;
df = @(x) x^2-1;

epsilon = 1.0e-8;%error

%initial
x1=[0.1 0.2 0.9 9.0];
x2=[-0.1,0.1; -0.2, 0.2; -2, 0.9; 0.9, 9.0];

%计算结果
fprintf("牛顿法\n");
%牛顿法-1
[y,m,list1] = Newton(f,df,x1(1),epsilon);
error = abs(list1 - 0);
fprintf("初值: x0=%f, 根: %.15f, 迭代次数: %d\n",x1(1),y,m);
fprintf("error[0]=%f\n",error(1));
for k = 1:length(list1)-1
    fprintf("error[%d]=%.15f; error[%d]/error[%d]^3=%.15f\n",k,error(k+1),k,k-1,error(k+1)/error(k)^3);
end
%牛顿法-2
[y,m,list1] = Newton(f,df,x1(2),epsilon);
error = abs(list1 - 0);
fprintf("初值: x0=%f, 根: %.15f, 迭代次数: %d\n",x1(2),y,m);
fprintf("error[0]=%f\n",error(1));
for k = 1:length(list1)-1
    fprintf("error[%d]=%.15f; error[%d]/error[%d]^3=%.15f\n",k,error(k+1),k,k-1,error(k+1)/error(k)^3);
end
%牛顿法-3
[y,m,list1] = Newton(f,df,x1(3),epsilon);
error = abs(list1 + sqrt(3));
fprintf("初值: x0=%f, 根: %.15f, 迭代次数: %d\n",x1(3),y,m);
fprintf("error[0]=%f\n",error(1));
for k = 1:length(list1)-1
    fprintf("error[%d]=%.15f; error[%d]/error[%d]^2=%.15f\n",k,error(k+1),k,k-1,error(k+1)/error(k)^2);
end
%牛顿法-4
[y,m,list1] = Newton(f,df,x1(4),epsilon);
error = abs(list1 - sqrt(3));
fprintf("初值: x0=%f, 根: %.15f, 迭代次数: %d\n",x1(4),y,m);
fprintf("error[0]=%f\n",error(1));
for k = 1:length(list1)-1
    fprintf("error[%d]=%.15f; error[%d]/error[%d]^2=%.15f\n",k,error(k+1),k,k-1,error(k+1)/error(k)^2);
end
fprintf("弦截法\n");
%弦截法-1
[z,m,list2] = Secant(f,x2(1,1),x2(1,2),epsilon);
fprintf("初值: x0=%f, x1=%f, 根: %.15f, 迭代次数: %d\n",x2(1,1),x2(1,2),z,m);
error = abs(list2 - 0);
fprintf("error[0]=%f\n",error(1));
for k = 1:length(list2)-1
    fprintf("error[%d]=%.15f\n",k,error(k+1));
end
%弦截法-2
[z,m,list2] = Secant(f,x2(2,1),x2(2,2),epsilon);
fprintf("初值: x0=%f, x1=%f, 根: %.15f, 迭代次数: %d\n",x2(2,1),x2(2,2),z,m);
error = abs(list2 - 0);
fprintf("error[0]=%f\n",error(1));
for k = 1:length(list2)-1
    fprintf("error[%d]=%.15f\n",k,error(k+1));
end
%弦截法-3
[z,m,list2] = Secant(f,x2(3,1),x2(3,2),epsilon);
fprintf("初值: x0=%f, x1=%f, 根: %.15f, 迭代次数: %d\n",x2(3,1),x2(3,2),z,m);
error = abs(list2 - sqrt(3));
fprintf("error[0]=%f\n",error(1));
for k = 1:length(list2)-1
    fprintf("error[%d]=%.15f; error[%d]/error[%d]^1.618=%.15f\n",k,error(k+1),k,k-1,error(k+1)/error(k)^1.618);
end
%弦截法-4
[z,m,list2] = Secant(f,x2(4,1),x2(4,2),epsilon);
fprintf("初值: x0=%f, x1=%f, 根: %.15f, 迭代次数: %d\n",x2(4,1),x2(4,2),z,m);
error = abs(list2 - sqrt(3));
fprintf("error[0]=%f\n",error(1));
for k = 1:length(list2)-1
    fprintf("error[%d]=%.15f; error[%d]/error[%d]^1.618=%.15f\n",k,error(k+1),k,k-1,error(k+1)/error(k)^1.618);
end


%通用函数

function [x,n,l] = Newton(f, df, x0 ,e)
n = 0;
l = [x0];
while 1
    n = n + 1;
    x = x0 - f(x0)/df(x0);
    l = [l x];
    if abs(x-x0)<e
        break;
    end
    x0 = x;
end
end


function [x,n,l] = Secant(f, x0, x1 ,e)
n = 0;
l = [x0 x1];
while 1
    n = n + 1;
    x = x1 - f(x1)/(f(x1)-f(x0))*(x1-x0);
    l = [l x];
    if abs(x-x1)<e
        break;
    end
    x0 = x1;
    x1 = x;
end
end
