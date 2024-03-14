f = @(x) sin(x)

Standard_value = integral(@(x)f(x),1,5)%精确值

T=ComTa(f,1,5);
S=ComSim(f,1,5);

%误差计算
T_error=abs(T-Standard_value);
S_error=abs(S-Standard_value);

%误差阶计算
T_order=zeros(1,12);
S_order=zeros(1,12);
for i=2:12
    T_order(i)=log(T_error(i-1)/T_error(i))/log(2);
    S_order(i)=log(S_error(i-1)/S_error(i))/log(2);
end

%打印
fprintf("Trapezoid\n");
for i=1:12
    fprintf("point_number=%d; Integral_value=%.8f; error=%.8f; order=%f\n",2^i+1,T(i),T_error(i),T_order(i));
end    
fprintf("Simpson\n");
for i=1:12
    fprintf("point_number=%d; Integral_value=%.16f; error=%.16f; order=%f\n",2^i+1,S(i),S_error(i),S_order(i));
end    


%复化梯形积分
function T=ComTa(f,a,b)
h=(b-a)./(2.^(1:12));
T=zeros(1,12);
T(1)=(b-a)*(f(a)+2*f((a+b)/2)+f(b))/4;
for j=2:12
  subtotal = 0;
  for i=1:2^(j-1)
    subtotal = subtotal + f(a+(2*i-1)*h(j));
  end
  T(j) = T(j-1)/2+h(j)*subtotal;
end
end


%复化Simpson积分公式
function S=ComSim(f,a,b)
h=(b-a)./(2.^(1:12));
S=zeros(1,12);
S(1)=(b-a)*(f(a)+4*f((a+b)/2)+f(b))/6;
for j=2:12
  S(j) = h(j)/3*(f(a)+f(b)+4*sum(f(a+h(j)*(2*(1:2^(j-1))-1)))+2*sum(f(a+h(j)*2*(1:(2^(j-1)-1)))));
end
end