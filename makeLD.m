function [L, D] = makeLD(Sigma)
    % `makeLD` Decompose a covariance matrix into the matrices `L` and `D`.
    %
    %   `[L, D] = makeLD(Sigma)` returns the inverse `L` of the lower Cholesky
    %   factor of `Sigma` and the diagonal matrix `D` holding the inverse of
    %   the diagonal of `L`, so that `D * L` has a unit diagonal.
    %
    %   ## Arguments
    %   - `Sigma` (matrix): Positive definite covariance matrix.
    %
    %   ## Returns
    %   - `L` (matrix): Inverse of the lower Cholesky factor of `Sigma`.
    %   - `D` (matrix): Diagonal matrix with entries `1 ./ diag(L)`.
    %
    %   ## Notes
    %   - This function is intended for internal use only.
    %
    %   See also `makeB`, `makeOmega`

    Linv = chol(Sigma, 'lower');
    L = inv(Linv);
    D = diag(1 ./ diag(L));
end
