A = [3, 2, 5, 4, 6 ;...
    2, 1, 3, -7, 8;...
    5, 3, 2, 5, -4;...
    4, -7, 5, 1, 3; ...
    6, 8,-4, 3, 8];
[D, V] = Jacobi(A,100, 1e-4);

n = length(A);
for k = 1:n
   fprintf("r%d = %.15f, ", k, D(k));
   fprintf("v%d = (%.15f",k, V(1,k));
   for j = 2:n
       fprintf(", %.15f",V(j,k));
   end
   fprintf(")\n");
end

function [D, V] = Jacobi(A, M, tol)
n = size(A, 1);
V = eye(n); %eigenvetors

for k = 1:M
    a_pq = 0;
    p = 1;
    q = 1;
    %find max A_pq
    for i = 1:n-1
        for j = i+1:n
            if abs(A(i,j))>a_pq
                a_pq = abs(A(i,j));
                p = i;
                q = j;
            end
        end
    end
    %compute tan
    if A(p,q)~=0
        s = (A(q,q) - A(p,p))/A(p,q)/2;
        if s >= 0
            t = 1/(s + sqrt(1+s^2));
        else
            t = 1/(s - sqrt(1+s^2));
        end
        c = 1/sqrt(1+t^2);
        d = c*t;
    else
        c = 1;
        d = 0;
    end
    
    %matrix Q
    Q = eye(n);
    Q(p,p) = c;
    Q(p,q) = d;
    Q(q,p) = -d;
    Q(q,q) = c;
    
    %rotate
    A = transpose(Q)*A*Q;
    V = V*Q;
    D = diag(A); %eigenvalues
    %test
    if norm(A-diag(D),"fro")<tol
        return;
    end
end
fprintf('Can not compute the eigenvalues of square matrix A within max iteration steps %d.\n', M);
end
