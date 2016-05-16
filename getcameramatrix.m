function P = getcameramatrix(x, X)

A = get_A(x, X);
[U S V] = svd(A);
P = V(:, size(V,2));
P = reshape(P, 4, 3)';

return