function val = logdet(Sigma)
    if ~isequal(Sigma, Sigma') || ~all(eig(Sigma) > 0)
        error("logdet: Matrix must be symmetric positive definite.")
    end
    R = chol(Sigma);
    val = 2 * sum(log(diag(R)));
end
