function [a, b] = for_back(f, lamda0, h0)
t = 2;
h = h0;
lamda1 = lamda0;
lamda = lamda1;
lamda2 = lamda1 + h;
% if f(lamda2) >= f(lamda1)
%     h = -h;
%     lamda1 = lamda2;
%     lamda2 = lamda1 + h;
% end
while f(lamda2) < f(lamda1)
    h = h*t;
    lamda = lamda1;
    lamda1 = lamda2;
    lamda2 = lamda1 + h;
end
a = min([lamda,lamda2]);
b = max([lamda,lamda2]);
end