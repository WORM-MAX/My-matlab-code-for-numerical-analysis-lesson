function x = guass_seidel(A, b, M, tol)
% Inputs: full or sparse matrix A,b
% number of iterations M, tolerance tol 
% Output: solution x
n = length(b);                          % find n
L = tril(A);                           % extract diagonal of A
R = A-tril(L);                        % R is the remainder
x=[0;0;0;0;0;0;0;0;0];                          % initialize vector x
for k=1:M                               % loop for iteration
    u = inv(L)*(b-R*x); 
    if norm(u-x, inf) < tol
        x = u;
        fprintf('k = %d, x(k) = (', k);
        for j=1:n-1
            fprintf('%.16f, ', u(j));
        end
        fprintf('%.16f), ', u(n));
        fprintf('||x(k) - x(k-1)|| = %.16f. \n', norm(u-x, inf));
        return;
    else
        x = u;
    end
end                                       % End of iteration loop

fprintf('can not solve the linear equation Ax = b.\n');
end