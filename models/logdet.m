function val = logdet(Sigma)
    % `logdet` Compute the log determinant of a covariance matrix.
    %
    %   `val = logdet(Sigma)` returns the log determinant of the
    %   symmetric positive definite matrix `Sigma`.
    %
    %   ## Arguments
    %   - `Sigma` (matrix): Symmetric positive definite matrix.
    %
    %   ## Returns
    %   - `val` (number): Log determinant of `Sigma`.
    if ~isequal(Sigma, Sigma') || ~all(eig(Sigma) > 0)
        error("logdet: Matrix must be symmetric positive definite.")
    end
    R = chol(Sigma);
    val = 2 * sum(log(diag(R)));
end
