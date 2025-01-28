function [L, D] = makeLD(Sigma)
    Linv = chol(Sigma, 'lower');
    L = inv(Linv);
    D = diag(1 ./ diag(L));
end
