function x = gauss_pivot(A, b)
% Gaussian elimination algorithm
% Inputs: A coefficients and b right side vector
% Output: the solution of Ax = b
assert(size(A, 1) == size(b, 1));
n = size(A, 1);

for k = 1:n-1
    %find largest pivot
    [~,maxix] = max(abs(A(k:n,k)));
    p = maxix + k -1;
    %Swap row
    A([k,p],:) = A([p,k],:);
    b([k,p]) = b([p,k]);
    for i = k+1:n
        z = A(i, k) / A(k, k);
        for j = k+1:n
            A(i, j) = A(i, j) - z * A(k, j);
        end
        b(i) = b(i) - z * b(k);
    end    
end

for i = n:-1:1
    tmp = 0;
    for j = i+1:n
        tmp = tmp + A(i, j) * x(j);
    end        
    x(i) = (b(i)  - tmp) / A(i, i);
end

end
