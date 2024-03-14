function [lam, x]=repowerit(A, x, M, tol)
% Power Iteration
% Computes dominant eigenvector of square matrix
% Input: matrix A, initial (nonzero) vector x, number of max iteration steps M
% Output: dominant eigenvalue lam, eigenvector x
assert(size(A, 1) == size(x, 1));
n = size(A, 1);
B = inv(A-1.2*eye(n));
v = x/norm(x); 
mu = v'*B*v;
for k = 1:M
    x = B*v;                                     % power step
    u = x/norm(x);                            % normalize vector
    lam = u'*B*u;                             % Rayleigh quotient approximation
    
    fprintf('k = %d, lam = %.16f, ', k, lam); 
    fprintf('u(k) = (', k);                     % output x(k)
    for j=1:n-1
        fprintf('%.16f, ', u(j));
    end
    fprintf('%.16f).\n', u(n));
    
    if abs(lam-mu) < tol || norm(u-v, 2) < tol
        x = u;
        return;
    else
        v = u;
        mu = lam;
    end
end
fprintf('Can not compute the dominant eigenvector of square matrix A within max iteration steps %d.\n', M);

end