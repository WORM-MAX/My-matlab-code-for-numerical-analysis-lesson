xi = [0.25; 0.5; 0.75; 1; 1.25; 1.5; 1.75; 2; 2.25; 2.5];
yi = [1.284; 1.648; 2.117; 2.718; 3.427; 2.798; 3.534; 4.456; 5.465; 5.894];

cof = linlsf(xi,yi);
fprintf("a = %f, b = %f, R^2 = %f",cof(1),cof(2),cof(3));
function cof = linlsf(xi,yi)
% Linear least square fitting
% Inputs: xi, yi sample points and values
% Output: coefficients of linear functions
assert(length(xi) == length(yi), 'the number points in xi and yi should be equal!');

A = [sin(xi) cos(xi)];
B = transpose(A);      
b = yi;
X = B*A \ (B*b);
cof = [X; transpose(b)*b-transpose(X)*B*A*X];

end