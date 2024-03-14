clear,clc

f=@(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
g=@(x)[-2*(1-x(1))-400*x(1)*(-x(1)^2+x(2)),200*(-x(1)^2+x(2))];%griend of f
H=@(x)[2+100*(8*x(1)^2-4*(-x(1)^2+x(2))),-400*x(1)%hessian matrix of f
    -400*x(1),200];
%fileID = fopen('data.txt','w');

fprintf('最速下降法:\n');
%fprintf(fileID,'最速下降法:\n');
gradient_descent_fminunc(f,g,[0 0],1E-4);

fprintf('\n\n');
%fprintf(fileID,'\n\n');

fprintf('牛顿法:\n');
%fprintf(fileID,'牛顿法:\n');
newton_fminunc(f,g,H,[0 0],1E-4);


function [x,fval]=gradient_descent_fminunc(f,g,x0,e,fileID)
%最速下降法求无约束非线性方程极小值
%输入 
%f:代求函数
%g:函数梯度
%x0:寻找x0附近的极小值点
%e:误差限度
%输出
%x: 极小值点
%fval: 极小值
    x=x0;g1=g(x);i=1;
    while(norm(g1)>e)
        f1=@(lambda)f(x-lambda*g1);%lambda 是单变量函数的自变量
        %这里虽然设置了lambda起始值是0，但输出的任然是正数，这是极小值原理保证的
        [lambda,~]=golden_section_fmin(f1,0);
        x=x-lambda*g1;
        g1=g(x);
        %fprintf('lamba=%f\n',lambda);
        fprintf('第%d次迭代：f(x_i)=%20.15e, x_1=%20.15e,x_2=%20.15e\n'...
            ,i,f(x),x(1),x(2));
        if(nargin==5)
            fprintf(fileID,'第%d次迭代：f(x_i)=%20.15e, x_1=%20.15e,x_2=%20.15e\n'...
                ,i,f(x),x(1),x(2));
        end
        i=i+1;
    end
    fval=f(x);
end

function [x,fval]=newton_fminunc(f,g,H,x0,e,fileID)
%牛顿法求无约束非线性方程极小值
%输入 
%f:待求函数
%g:函数梯度
%H:Hermitian matrix 
%x0:寻找x0附近的极小值点

%输出
%x: 极小值点
%fval: 极小值
    n=length(x0);%维度
    x=x0;g1=g(x);
    i=1;
    while(norm(g1)>e)
        p=-g1/H(x);
        f1=@(lambda)f(x+lambda*p);%lambda 是单变量函数的自变量
        [lambda,~]=golden_section_fmin(f1,0)
        x=x+lambda*p
        fprintf('第%d次迭代：f(x_i)=%20.15e, x_1=%20.15e,x_2=%20.15e\n'...
            ,i,f(x),x(1),x(2));
        if(nargin==6)
            fprintf(fileID,'第%d次迭代：f(x_i)=%20.15e, x_1=%20.15e,x_2=%20.15e\n'...
            ,i,f(x),x(1),x(2));
        end
        i=i+1;
        g1=g(x);
    end
    fval=f(x);
end

function [x,fval]=golden_section_fmin(f,x0)
%黄金分割法求极小值
%输入 
%f:代求函数
%x0:寻找x0附近的极小值点
%输出
%x: 极小值点
%fval: 极小值
    e0=1E-5;
    h=1;t=2;
    x2=x0+h;x1=x0;x=x0;
    %进退法找区间
    f(x2)
    f(x1)
    if(f(x2)>f(x1))
        h=-h;x1=x2;x2=x1+h;
    end
    while(f(x2)<f(x1))
        h=h*t;
        x=x1;x1=x2;x2=x1+h;
    end
    a=min(x,x2);b=max(x,x2);
    %黄金分割法迭代找极小值
    t=(-1+sqrt(5))/2;

    N=100000;%防止不收敛
    ak = b-t*(b-a);
    bk = a+t*(b-a);
    ya = f(ak);
    yb = f(bk);
    for k = 1:N
        if(ya <= yb)
            if(b-a<e0)
                x = ak;
                fval = ya;
                return;
            end
            %a = a;
            b = bk;
            bk = ak;%利用前一步算好的值
            ak = b-t*(b-a);
            yb = ya;%利用前一步算好的值
            ya = f(ak);
        else
            if(b-a<e0)
                x = bk;
                fval = yb;
                return;
            end
            %b = b;
            a = ak;
            ak = bk;%利用前一步算好的值     
            bk = a+t*(b-a);
            ya = yb;%利用前一步算好的值
            yb = f(bk);            
        end
    end
end